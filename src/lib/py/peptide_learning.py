""" Module for peptide learning """
from __future__ import annotations
from typing import TYPE_CHECKING
from pathlib import Path
import dataclasses

import pandas as pd
import numpy as np

import torch
import torch.nn as nn

from torch.utils.data import Dataset, DataLoader
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig


DEFAULT_SCALAR_FEATURES = [
    'N_Term_Cleaved',
    'C_Term_Cleaved',
    'Miscleavages',
    'Charge',
    'HyperScore',
    'weighted_spectral_entropy',
    'hypergeometric_probability',
    'delta_RT_loess',
    'intersection'
]

@dataclasses.dataclass
class Params:
    enzyme: str = 'trypsin'
    window_size: int = 16
    scaler_features: list[str] = dataclasses.field(default_factory=lambda: DEFAULT_SCALAR_FEATURES)


class PeptideHead(nn.Module):
    def __init__(self, embed_dim: int, hidden_dim: int = 512, n_scalar: int = 0):
        super().__init__()
        in_dim = 3 * embed_dim + n_scalar
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(p=0.2),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)  # logits

class ProteinPriorHead(nn.Module):
    def __init__(self, embed_dim: int, hidden_dim: int = 128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x.new_tensor(self.net(x).squeeze(-1))

class PeptideLearner:
    def __init__(self, unique_peptide_table: pd.DataFrame, shared_peptide_table: pd.DataFrame,
            database:dict[str, str], params: Params=None):
        self.unique_peptide_table = unique_peptide_table
        self.shared_peptide_table = shared_peptide_table
        self.database = database

        self.params = params or Params()
        self.head: PeptideHead | None = None
        self.esmc_device: torch.device | None = None
        self.esmc: ESMC | None = None
        self.esmc_embed_dim: int | None = None

        self._emb_cache: dict[str, torch.Tensor] = {}
        self.stage2_scalar_features = ["N_Term_Cleaved", "C_Term_Cleaved", "Miscleavages"]

    def _get_device(self) -> torch.device:
        if torch.backends.mps.is_available():
            return torch.device("mps")
        if torch.cuda.is_available():
            return torch.device("cuda")
        return torch.device("cpu")

    def init_esmc(self, model_name: str):
        device = self._get_device()
        self.esmc_device = device

        model = ESMC.from_pretrained(model_name)
        model = model.to(device)
        model.eval()
        for p in model.parameters():
            p.requires_grad = False

        self.esmc = model

        # Get embedding dim from model config (fallback to infer later)
        embed_dim = getattr(getattr(model, "config", None), "hidden_size", None)
        self.esmc_embed_dim = embed_dim

    @torch.no_grad()
    def _embed_sequence_esmc(self, seq: str) -> torch.Tensor:
        protein = ESMProtein(sequence=seq)
        protein_tensor = self.esmc.encode(protein).to(self.esmc_device)
        out = self.esmc.logits(protein_tensor, LogitsConfig(sequence=False, return_embeddings=True))

        # out.embeddings: (1, L, D)
        emb = out.embeddings.mean(dim=1).squeeze(0)  # (D,)
        return emb

    @torch.no_grad()
    def _get_scalar_features(self, row: pd.Series) -> torch.Tensor:
        scalars = torch.tensor(
            [
                row[f] for f in self.params.scaler_features
            ], dtype=torch.float32
        )
        return scalars

    @torch.no_grad()
    def build_feature_vector(self, row: pd.Series, include_scalars: bool = True) -> torch.Tensor:
        n_seq = self._normalize_sequence(row.get("N_Flanking"))
        pep = self._normalize_sequence(row.get("Peptide"))
        c_seq = self._normalize_sequence(row.get("C_Flanking"))

        e_n = self._embed_sequence_esmc(n_seq)
        e_p = self._embed_sequence_esmc(pep)
        e_c = self._embed_sequence_esmc(c_seq)

        if include_scalars:
            scalars = self._get_scalar_features(row)
            x = torch.cat([e_n, e_p, e_c, scalars], dim=0)
        else:
            x = torch.cat([e_n, e_p, e_c], dim=0)
        return x

    def _normalize_sequence(self, value: object | None) -> str:
        if value is None:
            return ''
        if isinstance(value, str):
            return value
        if pd.isna(value):
            return ''
        return str(value)

    def _cache_features_from_table(self, table: pd.DataFrame) -> dict[str, torch.Tensor]:
        assert table is not None, "A table must be provided to cache features."
        assert self.esmc is not None, "Call init_esmc() first."

        xs = []
        ys = []
        idxs = []
        scalars_list = []

        for idx, row in table.iterrows():
            y = float(row["Probability"])
            if not np.isfinite(y):
                continue
            y = max(0.0, min(1.0, y))  # clamp to [0,1]

            x = self.build_feature_vector(row, include_scalars=False)
            xs.append(x.cpu())
            ys.append(y)
            idxs.append(int(idx))

            scalars = self._get_scalar_features(row)
            scalars_list.append(scalars.cpu())

        X = torch.stack(xs, dim=0)  # (N, 3D)
        y = torch.tensor(ys, dtype=torch.float32)  # (N,)
        indices = torch.tensor(idxs, dtype=torch.long)  # (N,)
        scalars_cached = torch.stack(scalars_list, dim=0)  # (N, n_scalars)

        if getattr(self, "esmc_embed_dim", None) is None:
            self.esmc_embed_dim = X.shape[1] // 3

        cache = {
            "X": X,
            "y": y,
            "indices": indices,
            "scalars": scalars_cached,
            "embed_dim": self.esmc_embed_dim,
            "n_scalar": len(self.params.scaler_features)
        }

        return cache

    def cache_features(self, cache_path: Path | None = None):
        assert self.unique_peptide_table is not None

        cache = self._cache_features_from_table(self.unique_peptide_table)
        self.cached_X = cache["X"]
        self.cached_y = cache["y"]
        self.cached_indices = cache["indices"]
        self.cached_scalars = cache["scalars"]

        # If embed dim was not detected from config, infer from cached features
        if getattr(self, "esmc_embed_dim", None) is None:
            self.esmc_embed_dim = self.cached_X.shape[1] // 3

        if cache_path is not None:
            torch.save(cache, cache_path)

    def cache_features_for_table(self, table: pd.DataFrame, cache_path: Path | None = None) -> dict[str, torch.Tensor]:
        cache = self._cache_features_from_table(table)
        if cache_path is not None:
            torch.save(cache, cache_path)
        return cache

    def load_cached_features(self, cache_path: Path):
        obj = torch.load(cache_path, map_location="cpu")
        self.cached_X = obj["X"]
        self.cached_y = obj["y"]
        self.cached_indices = obj.get("indices")
        self.cached_scalars = obj.get("scalars")

        if getattr(self, "esmc_embed_dim", None) is None:
            self.esmc_embed_dim = self.cached_X.shape[1] // 3

    def train_head(
        self,
        epochs: int = 5,
        batch_size: int = 64,
        lr: float = 1e-3,
        weight_decay: float = 1e-4,
        hidden_dim: int = 512,
        val_fraction: float = 0.1,
        seed: int = 1231234
    ):
        assert hasattr(self, "cached_X") and hasattr(self, "cached_y"), "Run cache_features() first."

        device = self._get_device()

        # Shuffle/split
        torch.manual_seed(seed)
        n = self.cached_X.shape[0]
        idx = torch.randperm(n)
        n_val = int(round(val_fraction * n))
        val_idx = idx[:n_val]
        tr_idx = idx[n_val:]

        # Use cached scalars if available; otherwise compute from unique_peptide_table
        if hasattr(self, "cached_scalars") and self.cached_scalars is not None:
            scalars = self.cached_scalars
        else:
            if not hasattr(self, "cached_indices") or self.cached_indices is None:
                raise RuntimeError("Cached features missing indices; please re-run cache_features().")
            if self.unique_peptide_table is None:
                raise RuntimeError("unique_peptide_table required when scalars are not cached.")
            scalars_np = self.unique_peptide_table.loc[
                self.cached_indices.tolist(), self.params.scaler_features
            ].to_numpy(dtype=np.float32)
            scalars = torch.from_numpy(scalars_np)

        X_full = torch.cat([self.cached_X, scalars], dim=1)

        X_tr = X_full[tr_idx]
        y_tr = self.cached_y[tr_idx]
        X_val = X_full[val_idx]
        y_val = self.cached_y[val_idx]

        train_loader = DataLoader(
            list(zip(X_tr, y_tr)),
            batch_size=batch_size,
            shuffle=True
        )
        val_loader = DataLoader(
            list(zip(X_val, y_val)),
            batch_size=batch_size,
            shuffle=False
        )

        head = PeptideHead(
            embed_dim=self.esmc_embed_dim,
            hidden_dim=hidden_dim,
            n_scalar=len(self.params.scaler_features)
        ).to(device)
        loss_fn = nn.BCEWithLogitsLoss()
        opt = torch.optim.AdamW(head.parameters(), lr=lr, weight_decay=weight_decay)

        for epoch in range(1, epochs + 1):
            head.train()
            tr_losses = []

            for xb, yb in train_loader:
                xb = xb.to(device)
                yb = yb.to(device)

                opt.zero_grad()
                logits = head(xb)
                loss = loss_fn(logits, yb)
                loss.backward()
                opt.step()
                tr_losses.append(loss.item())

            head.eval()
            val_losses = []
            with torch.no_grad():
                for xb, yb in val_loader:
                    xb = xb.to(device)
                    yb = yb.to(device)
                    logits = head(xb)
                    loss = loss_fn(logits, yb)
                    val_losses.append(loss.item())

            print(
                f"Epoch {epoch:02d} | "
                f"train BCE={float(np.mean(tr_losses)):.4f} | "
                f"val BCE={float(np.mean(val_losses)):.4f}"
            )

        self.head = head

    @torch.no_grad()
    def predict_probability(self, row: pd.Series) -> float:
        assert self.head is not None, "Train the head first."
        device = next(self.head.parameters()).device
        x = self.build_feature_vector(row).unsqueeze(0).to(device)
        logit = self.head(x).squeeze(0)
        prob = torch.sigmoid(logit).item()
        return float(prob)

    @torch.no_grad()
    def predict_probabilities_batch(
        self,
        rows: pd.DataFrame,
        batch_size: int = 64,
        use_cache: bool = False,
        cache_path: Path | None = None,
        cache: dict[str, torch.Tensor] | None = None
    ) -> np.ndarray:
        """Predict probabilities for multiple peptides efficiently."""
        assert self.head is not None, "Train the head first."
        device = next(self.head.parameters()).device

        if cache is None and cache_path is not None and cache_path.exists():
            potential_cache = torch.load(cache_path, map_location="cpu")
            if (
                potential_cache.get("embed_dim") == self.esmc_embed_dim and
                potential_cache.get("n_scalar") == len(self.params.scaler_features)
            ):
                cache = potential_cache
            else:
                cache = None

        if cache is None and use_cache:
            cache = self.cache_features_for_table(rows, cache_path=cache_path)

        if cache is not None:
            X = cache["X"]
            scalars = cache["scalars"]
            X_full = torch.cat([X, scalars], dim=1)

            probs = []
            for i in range(0, X_full.shape[0], batch_size):
                xb = X_full[i:i + batch_size].to(device)
                logits = self.head(xb)
                batch_probs = torch.sigmoid(logits).cpu().numpy()
                probs.extend(batch_probs)

            return np.array(probs)

        probs = []
        for i in range(0, len(rows), batch_size):
            batch_rows = rows.iloc[i:i+batch_size]
            xs = torch.stack([
                self.build_feature_vector(row)
                for _, row in batch_rows.iterrows()
            ]).to(device)
            logits = self.head(xs)
            batch_probs = torch.sigmoid(logits).cpu().numpy()
            probs.extend(batch_probs)

        return np.array(probs)

    def apply_stage1_to_shared(self, batch_size: int = 64, cache_path: Path | None = None):
        """
        Apply the trained Stage 1 head to the shared peptide table
        and store the reliability weight per row.
        """
        assert self.head is not None, "Train Stage 1 head first."
        assert self.shared_peptide_table is not None

        probs = self.predict_probabilities_batch(
            self.shared_peptide_table,
            batch_size=batch_size,
            use_cache=True,
            cache_path=cache_path
        )

        self.shared_peptide_table = self.shared_peptide_table.copy()
        self.shared_peptide_table["Stage1Weight"] = probs

    def add_peptide_group_ids(self):
        """
        Assign an integer group id for each unique peptide.
        Required for grouped softmax in Stage 2.
        """
        assert self.shared_peptide_table is not None

        df = self.shared_peptide_table.copy()

        # Ensure peptide is normalized
        df["Peptide"] = df["Peptide"].astype(str)

        # Assign group ids
        peptide_to_gid = {
            pep: i
            for i, pep in enumerate(df["Peptide"].unique())
        }

        df["PeptideGroupID"] = df["Peptide"].map(peptide_to_gid)

        self.shared_peptide_table = df
        self.peptide_group_map = peptide_to_gid

    @torch.no_grad()
    def _embed_sequence_esmc(self, seq: str) -> torch.Tensor:
        seq = self._normalize_sequence(seq)

        if seq in self._emb_cache:
            return self._emb_cache[seq]

        protein = ESMProtein(sequence=seq)
        protein_tensor = self.esmc.encode(protein).to(self.esmc_device)
        out = self.esmc.logits(protein_tensor, LogitsConfig(sequence=False, return_embeddings=True))

        emb = out.embeddings.mean(dim=1).squeeze(0).detach().cpu()  # cache on CPU
        self._emb_cache[seq] = emb
        return emb

    @torch.no_grad()
    def build_stage2_feature_vector(self, row: pd.Series, include_scalars: bool = True) -> torch.Tensor:
        n_seq = self._normalize_sequence(row.get("N_Flanking"))
        pep = self._normalize_sequence(row.get("Peptide"))
        c_seq = self._normalize_sequence(row.get("C_Flanking"))

        e_n = self._embed_sequence_esmc(n_seq)
        e_p = self._embed_sequence_esmc(pep)
        e_c = self._embed_sequence_esmc(c_seq)

        x = torch.cat([e_n, e_p, e_c], dim=0)

        if include_scalars:
            scalars = torch.tensor(
                [float(row.get(f, 0.0)) for f in self.stage2_scalar_features],
                dtype=torch.float32
            )
            x = torch.cat([x, scalars], dim=0)

        return x

    def prepare_stage2_training_tensors(self) -> None:
        """
        Build and store tensors needed for Stage-2 EM training:
        - stage2_X: (N, D) feature matrix for shared_peptide_table rows
        - stage2_group_ids: (N,) peptide group ids
        - stage2_weights: (N,) Stage-1 reliability weights
        """
        assert self.shared_peptide_table is not None, "shared_peptide_table is required."
        assert "PeptideGroupID" in self.shared_peptide_table.columns, "Run add_peptide_group_ids() first."
        assert "Stage1Weight" in self.shared_peptide_table.columns, "Run apply_stage1_to_shared() first."
        assert self.esmc is not None, "Call init_esmc() first."

        df = self.shared_peptide_table

        xs = []
        group_ids = []
        weights = []

        for _, row in df.iterrows():
            w = float(row["Stage1Weight"])
            if not np.isfinite(w):
                continue

            x = self.build_stage2_feature_vector(row, include_scalars=True)
            xs.append(x)

            group_ids.append(int(row["PeptideGroupID"]))
            weights.append(w)

        X = torch.stack(xs, dim=0)  # CPU tensors
        group_ids_t = torch.tensor(group_ids, dtype=torch.long)
        weights_t = torch.tensor(weights, dtype=torch.float32).clamp(0.0, 1.0)

        self.stage2_X = X
        self.stage2_group_ids = group_ids_t
        self.stage2_weights = weights_t

        embed_dim = X.shape[1] - len(self.stage2_scalar_features)
        self.stage2_embed_dim = embed_dim // 3
        self.stage2_n_scalar = len(self.stage2_scalar_features)

        print(
            f"Stage2 tensors built: N={X.shape[0]}, "
            f"embed_dim={self.stage2_embed_dim}, n_scalar={self.stage2_n_scalar}, "
            f"n_groups={int(group_ids_t.unique().numel())}"
        )

    def init_stage2_responsibilities(self):
        """
        Initialize EM responsibilities r_{p,i} uniformly within each peptide group.
        """
        assert hasattr(self, "stage2_group_ids"), "Run prepare_stage2_training_tensors() first."

        group_ids = self.stage2_group_ids
        N = group_ids.shape[0]

        r = torch.zeros(N, dtype=torch.float32)

        for g in torch.unique(group_ids):
            mask = group_ids == g
            k = mask.sum().item()
            r[mask] = 1.0 / k

        self.stage2_r = r

    def init_stage2_head(self, hidden_dim: int = 512):
        device = self._get_device()

        self.stage2_head = PeptideHead(
            embed_dim=self.stage2_embed_dim,
            hidden_dim=hidden_dim,
            n_scalar=self.stage2_n_scalar
        ).to(device)

    def stage2_m_step(
        self,
        epochs: int = 3,
        lr: float = 1e-3,
        weight_decay: float = 1e-4,
        batch_size: int = 1024
    ):
        device = self._get_device()

        X = self.stage2_X.to(device)
        group_ids = self.stage2_group_ids.to(device)
        r = self.stage2_r.to(device)
        w = self.stage2_weights.to(device)

        opt = torch.optim.AdamW(
            [
                {"params": self.stage2_head.parameters(), "lr": 1e-3},
                {"params": self.stage2_prior_head.parameters(), "lr": 1e-4},
            ],
            weight_decay=weight_decay
        )

        for epoch in range(1, epochs + 1):
            self.stage2_head.train()
            opt.zero_grad()

            context_logits = self.stage2_head(X)
            prior_logits = self.stage2_prior_head(self.stage2_protein_X.to(device))

            logits = context_logits + prior_logits
            probs = grouped_softmax(logits, group_ids)

            # EM weighted cross-entropy
            # Numerical safety constant to avoid zero.
            eps = 1e-8
            loss = -torch.sum(w * r * torch.log(probs + eps)) / w.sum()

            loss.backward()
            opt.step()

            print(f"[Stage2 M-step] epoch {epoch} | loss = {loss.item():.4f}")

    @torch.no_grad()
    def stage2_e_step(self):
        """
        E-step: update responsibilities r_{p,i} using the current Stage-2 head.
        """
        device = self._get_device()

        X = self.stage2_X.to(device)
        group_ids = self.stage2_group_ids.to(device)

        self.stage2_head.eval()

        context_logits = self.stage2_head(X)
        prior_logits   = self.stage2_prior_head(self.stage2_protein_X.to(device))

        logits = context_logits + prior_logits
        probs  = grouped_softmax(logits, group_ids)

        # Store updated responsibilities on CPU
        self.stage2_r = probs.detach().cpu()

    def cache_protein_embeddings(self, cache_path: Path | None = None):
        """
        Compute and cache ESM embeddings for full protein sequences.
        Used ONLY for protein-level priors. Optionally serialize the cache.
        """
        assert self.esmc is not None
        assert self.database is not None

        self.protein_emb_cache = {}

        for prot_id, seq in self.database.items():
            emb = self._embed_sequence_esmc(seq)  # uses cache internally
            self.protein_emb_cache[prot_id] = emb.cpu()

        print(f"Cached protein embeddings: {len(self.protein_emb_cache)}")

        if cache_path is not None:
            torch.save(self.protein_emb_cache, cache_path)

    def load_protein_embeddings(self, cache_path: Path):
        """
        Load cached protein embeddings from disk.
        """
        self.protein_emb_cache = torch.load(cache_path, map_location="cpu")
        print(f"Loaded protein embeddings: {len(self.protein_emb_cache)}")

    def prepare_stage2_protein_prior_tensors(self):
        """
        Build protein-embedding tensor aligned with stage2_X rows.
        """
        assert hasattr(self, "stage2_group_ids")
        assert hasattr(self, "protein_emb_cache")

        prot_embs = []

        for _, row in self.shared_peptide_table.iterrows():
            prot_id = row["ProteinAccession"]
            prot_embs.append(self.protein_emb_cache[prot_id])

        self.stage2_protein_X = torch.stack(prot_embs, dim=0)

    def init_stage2_prior_head(self):
        embed_dim = next(iter(self.protein_emb_cache.values())).shape[0]
        device = self._get_device()

        self.stage2_prior_head = ProteinPriorHead(
            embed_dim=embed_dim,
            hidden_dim=128
        ).to(device)

def grouped_softmax(logits: torch.Tensor, group_ids: torch.Tensor) -> torch.Tensor:
    out = torch.empty_like(logits)
    for g in torch.unique(group_ids):
        mask = group_ids == g
        out[mask] = torch.softmax(logits[mask], dim=0)
    return out


def entropy_per_group(group_ids, r, normalize=True, base=2, eps=1e-12):
    """
    group_ids: 1D array-like of ints (length N)
    r: 1D array-like of floats (length N), responsibilities per row
    """
    group_ids = np.asarray(group_ids)
    r = np.asarray(r, dtype=float)

    # clamp for numerical stability
    r = np.clip(r, eps, 1.0)

    df = pd.DataFrame({"group": group_ids, "r": r})

    def _entropy(x):
        if base == 2:
            h = -np.sum(x * np.log2(x))
        else:
            h = -np.sum(x * np.log(x)) / np.log(base)
        if normalize:
            k = len(x)
            if k <= 1:
                return 0.0
            h_max = np.log2(k) if base == 2 else np.log(k) / np.log(base)
            return float(h / h_max)
        return float(h)

    out = (
        df.groupby("group")["r"]
          .apply(lambda s: _entropy(s.to_numpy()))
          .reset_index(name="entropy")
    )

    # Also include group size for interpretation
    out["k"] = df.groupby("group").size().values
    return out

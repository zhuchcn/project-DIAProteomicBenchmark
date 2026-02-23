import AddContaminants
import CombineFasta
import ComputeFDR
import Convert2MzML
import DownloadData
from EcoliCustomDatabase.command import EcoliCustomDatabase
import FilterCanonicalProteins
import PhilosopherDatabase
import PhilosopherFilter
import PreparePeptideTrainingData
import RunFragPipe

SUBCOMMANDS = [
    AddContaminants.AddContaminants,
    CombineFasta.CombineFasta,
    ComputeFDR.ComputeFDR,
    Convert2MzML.Convert2MzML,
    DownloadData.DownloadData,
    EcoliCustomDatabase,
    FilterCanonicalProteins.FilterCanonicalProteins,
    PhilosopherDatabase.PhilosopherDatabase,
    PhilosopherFilter.PhilosopherFilter,
    PreparePeptideTrainingData.PreparePeptideTrainingData,
    RunFragPipe.RunFragPipe,
]

#!/usr/bin/env python3
from csv import DictReader
from os.path import isfile
from datetime import date
from dataclasses import dataclass, field
from typing import Tuple

@dataclass
class WaterYearSample:
    sampling_site: str
    replicate: str
    sampling_week: str
    depth: Tuple[int]
    depth_str: str
    date: str
    metadata_file_path: str
    core_GWC: float = field(default=None, init=False)
    core_pH: float = field(default=None, init=False)

    def __post_init__(self):
        self.date = date.fromisoformat(self.date)
        self.sample_id = self.__repr__()
        self.oneword_id = self.sample_id.replace('-','').replace('_','')
        if not isfile(self.metadata_file_path):
            raise FileNotFoundError(f"Metadata file {self.metadata_file_path} not found!")
        if self.core_GWC is None:
            self.core_GWC = self._get_GSC_property(self.metadata_file_path, 'GWC %')
        if self.core_pH is None:
            self.core_pH = self._get_GSC_property(self.metadata_file_path, 'pH')

    def __repr__(self):
        return f"{self.sampling_site}{self.replicate}{self.sampling_week}_{self.depth_str}"

    def _get_GSC_property(self, filePath, property):
        '''Get a property of a sample from the metadata file'''
        with open(filePath) as csvfile:
            reader = DictReader(csvfile)
            for row in reader:
                if row['sample_id'] == f"{self.sample_id}_GSC":
                    if property not in row.keys():
                        raise KeyError(f"Property {property} not found in metadata file at row {row['sample_id']}!")
                    if row[property] == '':
                        raise ValueError(f"Property {property} is empty in metadata file at row {row['sample_id']}!")
                    return float(row[property])

@dataclass
class UnfracSIPSample(WaterYearSample):
    isotope: str
    incubation_length: int
    grams_soil_extracted: float
    concentration_DNA_extracted: float
    total_ul_DNA_extracted: float
    h20_added: float

    def __post_init__(self):
        self.sample_id = self.__repr__()
        self.oneword_id = self.sample_id.replace('-','').replace('_','')

    def __repr__(self):
        return f"{self.sampling_site}{self.replicate}{self.sampling_week}_{self.depth[0]}-{self.depth[1]}_{self.isotope}O-{self.incubation_length}"
    

@dataclass
class SIPSample(UnfracSIPSample):
    isotope: str
    incubation_length: int
    grams_soil_extracted: float
    concentration_DNA_extracted: float
    total_ul_DNA_extracted: float
    ul_DNA_SIP_loaded: float
    h20_added: float
    plate: str
    tube: str
    dna_yield: float
    concentrations: "list[float]" = field(default_factory=list, repr=False)
    densities: "list[float]" = field(default_factory=list, repr=False)
    wells: "list[str]" = field(default_factory=list, repr=False)

    def __post_init__(self):
        self.sample_id = self.__repr__()
        self.oneword_id = self.sample_id.replace('-','').replace('_','')

    def __repr__(self):
        return f"{self.sampling_site}{self.replicate}{self.sampling_week}_{self.depth[0]}-{self.depth[1]}_{self.isotope}O-{self.incubation_length}"



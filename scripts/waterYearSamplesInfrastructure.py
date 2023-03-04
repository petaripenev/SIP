#!/usr/bin/env python3
from csv import DictReader
from os.path import isfile
from datetime import date
from dataclasses import dataclass, field, fields
from typing import Tuple
from warnings import warn

@dataclass
class WaterYearSample:
    id: int
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
        self.core_id = f"{self.sampling_site}{self.replicate}{self.sampling_week}_{self.depth_str}"
        self.oneword_id = self.core_id.replace('-','').replace('_','')
        if not isfile(self.metadata_file_path):
            raise FileNotFoundError(f"Metadata file {self.metadata_file_path} not found!")
        if self.core_GWC is None:
            self.core_GWC = self._get_experiment_property(self.metadata_file_path, self.core_id, 'GSC','GWC %')
        if self.core_pH is None:
            self.core_pH = self._get_experiment_property(self.metadata_file_path, self.core_id, 'GSC', 'pH')

    def __repr__(self):
        return f"{self.sampling_site}{self.replicate}{self.sampling_week}_{self.depth_str}"

    def _get_experiment_property(self, filePath, core_id, experiment, property):
        '''Get a property of a sample from the metadata file, 
        given the experiment name and the property name.
        Assumes the property is a float.'''
        with open(filePath) as csvfile:
            reader = DictReader(csvfile)
            for row in reader:
                if row['sample_id'] != f"{core_id}_{experiment}":
                    continue
                if property not in row.keys():
                    raise KeyError(f"Property {property} not found in metadata file at row {row['sample_id']}!")
                if row[property] == '':
                    warn(f"Property {property} is empty in metadata file at row {row['sample_id']}!", RuntimeWarning)
                    return None
                return float(row[property])
        raise KeyError(f"Sample {core_id}_{experiment} not found in metadata file!")

@dataclass
class DNAextractionSample(WaterYearSample):
    experiment: str
    grams_soil_extracted: float = field(default=None, init=False, metadata={'property': 'Soil for DNA extraction (gr)'})
    concentration_DNA_extracted: float = field(default=None, init=False, metadata={'property': 'DNA concentration (ng/ul)'})
    total_ul_DNA_extracted: float = field(default=None, init=False, metadata={'property': 'Volume DNA extraction (ul)'})

    def __post_init__(self):
        super().__post_init__()
        none_attributes = [field for field in fields(self) if getattr(self, field.name) is None]
        for attribute in none_attributes:
            if 'property' not in attribute.metadata.keys():
                warn(f"Attribute {attribute.name} is empty in metadata file for core {self.core_id}!", RuntimeWarning)
                continue
            property = attribute.metadata['property']
            setattr(self, attribute.name, self._get_experiment_property(self.metadata_file_path, self.core_id, self.experiment, property))

    def __repr__(self):
        return f"{self.sampling_site}{self.replicate}{self.sampling_week}_{self.depth_str}_{self.experiment}"

@dataclass
class FractionatedSIPSample(DNAextractionSample):
    ul_DNA_SIP_loaded: float
    plate: str
    tube: str
    dna_yield: float
    H2O_added: float = field(default=None, init=False, metadata={'property': 'H2O amount (ml)'})
    concentrations: "list[float]" = field(default_factory=list, repr=False)
    densities: "list[float]" = field(default_factory=list, repr=False)
    wells: "list[str]" = field(default_factory=list, repr=False)
    fraction_volumes: "list[float]" = field(default_factory=list, repr=False)

    def __post_init__(self):
        super().__post_init__()
        self.sample_id = self.__repr__()
        self.oneword_id = self.sample_id.replace('-','').replace('_','')
        self.area = [a*b for a,b in zip(self.densities,self.concentrations)]
        self.weighted_mean_density = sum(self.area)/sum(self.concentrations)
        self.remainingDNAperFraction = [((abs(conc)+conc)/2)*vol for conc,vol in zip(self.concentrations, self.fraction_volumes)]

    def __repr__(self):
        return f"{self.sampling_site}{self.replicate}{self.sampling_week}_{self.depth_str}_{self.experiment}"



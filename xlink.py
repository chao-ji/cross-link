import Peptide from peptide
import itertools

class XLink:
	def __init__(self, peptide1, peptide2, positions, charge, mass, param):
		self.pep_pair = [peptide1, peptide2]
		self.positions = positions
		self.charge = charge
		self.linker_mass = param['linker_mass']
		self.mol_weight = self.pep_pair[0].precursor_mass + self.pep_pair[1].precusor_mass + self.linker_mass
		
	def adjust_mass_list_per_peptide(self, index_this, index_that, mass, param):
		mass_list = self.pep_pair[index_this].get_mass_list(mass, param)
		link_pos = self.positions[index_this]
		precursor_mass = self.pep_pair[index_that].precursor_mass
		linker_mass = self.linker_mass
		length = self.pep_pair[index_this].length
		use_a_ion = param['use_a_ion']

		for i in range(length - 1):
			if i < link_pos:
				mass_list[i]['y'] += (linker_mass + precursor_mass)
			else:
				mass_list[i]['b'] += (linker_mass + precursor_mass)
				if use_a_ion:
					mass_list[i]['a'] += (linker_mass + precursor_mass)

		return massList	

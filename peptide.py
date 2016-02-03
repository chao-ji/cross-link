class Peptide:
	def __init__(self, sequence, protein_id, modification, mass, is_nterm):
		self.sequence = sequence
		self.length = len(self.sequence)
		self.protein_id = [protein_id]
		self.modification = modification
		self.mass_array = self.get_mass_array(mass)
		self.total_residue_mass = self.get_total_residue_mass()
		self.precursor_mass = self.get_precursor_mass(mass)
		self.is_nterm = is_nterm

	def get_mass_array(self, mass):
		mass_array = []
		sequence = self.sequence
		position = self.modification['position']
		delta_mass = self.modification['delta_mass']
		
		for i in range(self.length):
			mass_array.append(mass[sequence[i].upper()])
			
		for i in range(len(position)):
			mass_array[position[i]] += delta_mass[i]
			
		return mass_array
		
	def get_total_residue_mass(self):
		total_residue_mass = sum(self.mass_array)
		return total_residue_mass
		
	def get_precursor_mass(self, mass):
		precursor_mass = self.total_residue_mass
		precursor_mass = precursor_mass + mass['Hatom'] * 2 + mass['Oatom']
		return pm
		
	def get_mass_list(self, mass, param):
		mass_list = []
		fwd_mass = 0
		mass_array = self.mass_array
		total_residue_mass = self.total_residue_mass
		use_a_ion = param['use_a_ion']

		for i in range(self.length - 1):
			fragment = dict()
			fwd_mass += mass_array[i]
			rev_mass = total_residue_mass - fwd_mass
			fragment['b'] = fwd_mass + mass['b_ion_res']
			fragment['y'] = rev_mass + mass['y_ion_res']
			
			if use_a_ion:
				fragment['a'] = fwd_mass + mass['a_ion_res']
				
			mass_list.append(fragment)
			
		return mass_list	

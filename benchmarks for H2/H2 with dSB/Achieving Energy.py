#create molecule
name = 'H2'
charge = 0
multiplicity = 1
basis = 'sto-3g'
energy1 = []
energy2 = []
energy3 = []
energy4 = []
energy5 = []
hf_energies = []
exact_energies = []
Tx = 0
bond_lengths = np.linspace(0.2, 4.0, 50)
for bond_length in bond_lengths:
  geometry = get_molGeometry(name, bond_length)
  molecule = MolecularData(
    geometry=geometry,
    basis=basis,
    multiplicity=multiplicity,
    charge=charge)
  #run RHF calculations
  molecule = run_pyscf(molecule, run_scf=True)
  hf_energy = float(molecule.hf_energy)
  hf_energies.append(hf_energy)
  hf_data = molecule._pyscf_data['scf']
  #define active space
  n_active_electrons = 2
  n_active_orbitals = 2
  occupied_indices, active_indices = get_active_space(molecule, n_active_electrons, n_active_orbitals)
  #convert to fermionic Hamiltonian
  molecular_H = molecule.get_molecular_hamiltonian(occupied_indices=occupied_indices, active_indices=active_indices)
  if molecular_H[()] == None:
      molecular_H[()] = 0
  fermionic_H = get_fermion_operator(molecular_H)
  #add penalty term to ensure correct number of electrons in ground state
  weight = 5
  penalty_term = FermionOperator('', n_active_electrons)

  for i in range(molecular_H.n_qubits):
     penalty_term += FermionOperator(str(i)+'^ '+str(i), -1)
  fermionic_H += weight*penalty_term**2

  #convert to Pauli operator Hamiltonian
  binary_code = bravyi_kitaev_code(molecular_H.n_qubits)
  qubit_H = binary_code_transform(fermionic_H, binary_code)
  qubit_H.compress()

  #apply symmetry reductions and calculate minimum eigenvalue
  sectors = taper_qubits(qubit_H)
  qubit_H, min_eigenvalue = sector_with_ground(sectors)
  #consider exact energies as CASCI
  exact_energies.append(min_eigenvalue)
  #exact dissociatation limit calculation
  _, E_diss = sector_with_ground(sectors)
  m = count_qubits(qubit_H)

  ### XBK method ###

  #set r value
  r = 2

  #construct qubit Hamiltonians and C terms for XBK method
  qubit_Hs, qubit_Cs = [],[]
  for p in range(int(math.ceil(r/2+1))):
      qubit_Hs += [XBK_transform(qubit_H, r, p)]
      qubit_Cs += [construct_C(m, r, p)]

  #run XBK method
  XBK_energy, ground_state, t = XBK_dSB(qubit_Hs, qubit_Cs, r, starting_lam=0, num_steps=50, strength=1e3, verbose=False)
  Tx += t
  energy2.append(XBK_energy)
  #set r value
  r = 3

  #construct qubit Hamiltonians and C terms for XBK method
  qubit_Hs, qubit_Cs = [],[]
  for p in range(int(math.ceil(r/2+1))):
      qubit_Hs += [XBK_transform(qubit_H, r, p)]
      qubit_Cs += [construct_C(m, r, p)]

  #run XBK method
  XBK_energy, ground_state, t = XBK_dSB(qubit_Hs, qubit_Cs, r, starting_lam=0, num_steps=50, strength=1e3, verbose=False)
  Tx += t
  energy3.append(XBK_energy)

  #set r value
  r = 4

  #construct qubit Hamiltonians and C terms for XBK method
  qubit_Hs, qubit_Cs = [],[]
  for p in range(int(math.ceil(r/2+1))):
      qubit_Hs += [XBK_transform(qubit_H, r, p)]
      qubit_Cs += [construct_C(m, r, p)]

  #run XBK method
  XBK_energy, ground_state, t = XBK_dSB(qubit_Hs, qubit_Cs, r, starting_lam=0, num_steps=50, strength=1e3, verbose=False)
  Tx += t
  energy4.append(XBK_energy)

  #set r value
  r = 5

  #construct qubit Hamiltonians and C terms for XBK method
  qubit_Hs, qubit_Cs = [],[]
  for p in range(int(math.ceil(r/2+1))):
      qubit_Hs += [XBK_transform(qubit_H, r, p)]
      qubit_Cs += [construct_C(m, r, p)]

  #run XBK method
  XBK_energy, ground_state, t = XBK_dSB(qubit_Hs, qubit_Cs, r, starting_lam=0, num_steps=50, strength=1e3, verbose=False)
  Tx += t
  energy5.append(XBK_energy)

  #set r value
  r = 6

  #construct qubit Hamiltonians and C terms for XBK method
  qubit_Hs, qubit_Cs = [],[]
  for p in range(int(math.ceil(r/2+1))):
      qubit_Hs += [XBK_transform(qubit_H, r, p)]
      qubit_Cs += [construct_C(m, r, p)]

  #run XBK method
  XBK_energy, ground_state, t = XBK_dSB(qubit_Hs, qubit_Cs, r, starting_lam=0, num_steps=50, strength=1e3, verbose=False)
  Tx += t
  energy1.append(XBK_energy)

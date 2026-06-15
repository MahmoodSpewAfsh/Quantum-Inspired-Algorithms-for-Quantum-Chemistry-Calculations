# Create molecule
name = 'H2'
charge = 0
multiplicity = 1
basis = 'sto-3g'

r_values = [2, 3, 4, 5, 6]
all_energies = {r: [] for r in r_values}
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

    # Run RHF calculations
    molecule = run_pyscf(molecule, run_scf=True)
    hf_energies.append(float(molecule.hf_energy))

    # Define active space
    n_active_electrons = 2
    n_active_orbitals = 2
    occupied_indices, active_indices = get_active_space(molecule, n_active_electrons, n_active_orbitals)

    # Convert to fermionic Hamiltonian
    molecular_H = molecule.get_molecular_hamiltonian(
        occupied_indices=occupied_indices,
        active_indices=active_indices)
    if molecular_H[()] == None:
        molecular_H[()] = 0
    fermionic_H = get_fermion_operator(molecular_H)

    # Add penalty term
    weight = 5
    penalty_term = FermionOperator('', n_active_electrons)
    for i in range(molecular_H.n_qubits):
        penalty_term += FermionOperator(str(i)+'^ '+str(i), -1)
    fermionic_H += weight * penalty_term**2

    # Convert to Pauli operator Hamiltonian
    binary_code = bravyi_kitaev_code(molecular_H.n_qubits)
    qubit_H = binary_code_transform(fermionic_H, binary_code)
    qubit_H.compress()

    # Apply symmetry reductions and calculate minimum eigenvalue
    sectors = taper_qubits(qubit_H)
    qubit_H, min_eigenvalue = sector_with_ground(sectors)
    exact_energies.append(min_eigenvalue)

    # Exact dissociation limit (final value used in plot)
    _, E_diss = sector_with_ground(sectors)
    m = count_qubits(qubit_H)

    # XBK method — loop over r values
    for r in r_values:
        qubit_Hs, qubit_Cs = [], []
        for p in range(int(math.ceil(r / 2 + 1))):
            qubit_Hs.append(XBK_transform(qubit_H, r, p))
            qubit_Cs.append(construct_C(m, r, p))

        XBK_energy, ground_state, t = XBK_dSB(
            qubit_Hs, qubit_Cs, r,
            starting_lam=0, num_steps=50, strength=1e3, verbose=False)
        Tx += t
        all_energies[r].append(XBK_energy)

# Map r values to the energy lists
energy2, energy3, energy4, energy5, energy1 = [all_energies[r] for r in r_values]

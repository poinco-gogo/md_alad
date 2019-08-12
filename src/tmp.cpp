	// set 1-4 pair lists
	for (int i = 0; i < DihedralVector.size(); i++)
	{
		Atom* atom1 = DihedralVector[i].ptr_atom1;
//		Atom* atom2 = DihedralVector[i].ptr_atom2;
//		Atom* atom3 = DihedralVector[i].ptr_atom3;
		Atom* atom4 = DihedralVector[i].ptr_atom4;
/*
		if (atom1 -> PSFIndex == 5)// || atom4 -> PSFIndex == 5)
		{
			cout << "DEBUG hit at Dihedral i = " << i << endl;
			cout << 
			"periodicity = " << DihedralVector[i].n << endl;
//			atom1->PSFIndex << " " <<
//			atom2->PSFIndex << " " <<
//			atom3->PSFIndex << " " <<
//			atom4->PSFIndex << " " << endl;
		}*/

		atom1 -> scaled1_4Vector.push_back(atom4 -> PSFIndex);
		atom4 -> scaled1_4Vector.push_back(atom1 -> PSFIndex);
	}

bool PSF::set_LJ_parm(vector<Atom>& LJParmVector)
{
	if (!LJParmVector.size())
	{
		cerr << "error: please set LJ parm\n";
		return false;
	}
	
	for (int i = 0; i < ptr_atomVector -> size(); i++)
	{
		Atom& atom = ptr_atomVector -> at(i);
		for (int j = 0; j < LJParmVector.size(); j++)
		{
			Atom& LJatom = LJParmVector.at(j);
			//cout << atom->PSFAtomName << "  " 
			//	<< LJatom->PSFAtomName << '\n';
			if (atom.PSFAtomName == LJatom.PSFAtomName)
			{
				atom.epsilon = LJatom.epsilon;
				atom.Rmin_div2 = LJatom.Rmin_div2;
				atom.eps1_4 = LJatom.eps1_4;
				atom.Rmin1_4 = LJatom.Rmin1_4;
			
				break;
			}
			if (j == LJParmVector.size() - 1)
			{
				cerr << "error: LJ parameter missing for: "
					<< atom.PSFAtomName << '\n';
				return false;
			}
		}
	}

	// remove all elements
	LJParmVector.clear();
	cout << "REMARK LJ parameters are successfully assigned.\n";
	return true;

} // end of function() set_LJ_parm

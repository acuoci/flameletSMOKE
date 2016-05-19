{
	scalar m_in, m_out;
	scalar v_in, v_out;
	scalar H_in, H_out;

	m_in  = 0.;
	v_in  = 0.;
	H_in  = 0.;
	m_out = 0.;
	v_out = 0.;
	H_out = 0.;

	outInfo << runTime.value() << "\t";

	forAll(mesh.boundary(), patchi)
	{

		forAll(inletPatches, iP)
		{
			if(mesh.boundary()[patchi].name() == inletPatches[iP])
			{
				scalar flow =
				sum
				(
					(
						rho.boundaryField()[patchi]
						* U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
					)
				);

			   scalar Hflow =
				sum
				(
					(
						rho.boundaryField()[patchi]
						* H.boundaryField()[patchi]
						* U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
					)
				);

				scalar volFlow =
				sum
				(
					(
						U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
					)
				);

				m_in += mag(flow);
				v_in += mag(volFlow);
				H_in += mag(Hflow);
			}
		}

		forAll(outletPatches, oP)
		{
			if(mesh.boundary()[patchi].name() == outletPatches[oP])
			{
				scalar flow =
				sum
				(
					(
						rho.boundaryField()[patchi]
						* U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
					)
				);

				scalar Hflow =
				sum
				(
					(
						rho.boundaryField()[patchi] * H.boundaryField()[patchi]
						* U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
					)
				);

				scalar volFlow =
				sum
				(
					(
						U.boundaryField()[patchi]
						& mesh.boundary()[patchi].Sf()
					)
				);

				m_out += mag(flow);
				v_out += mag(volFlow);
				H_out += mag(Hflow);
			}
		}
	}

	outInfo << m_in << "\t" <<  v_in << "\t" <<  H_in << "\t";
	outInfo << m_out << "\t" <<  v_out << "\t" <<  H_out << "\t";
	outInfo << m_out-m_in << "\t" <<  v_out-v_in << "\t" <<  H_out-H_in << "\t";
	outInfo << endl;
}

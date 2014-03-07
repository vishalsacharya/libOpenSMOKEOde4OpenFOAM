
OpenSMOKEOde::OpenSMOKEOde(const Foam::ODESystem& openFOAM_odeSystem) :
openFOAM_odeSystem_(openFOAM_odeSystem)
{
	number_of_equations_ = openFOAM_odeSystem_.nEqns();
	yFoam_.resize(number_of_equations_);
	dyFoam_.resize(number_of_equations_);
}

int OpenSMOKEOde::Equations(const double t, const Eigen::VectorXd& y, Eigen::VectorXd& dy)
{	
	for(unsigned int i=0;i<number_of_equations_;i++)
		yFoam_[i] = y(i);
	openFOAM_odeSystem_.derivatives(t,yFoam_,dyFoam_);
	for(unsigned int i=0;i<number_of_equations_;i++)
		dy(i) = dyFoam_[i];

	return 0;
}

int OpenSMOKEOde::Print(const double t, const Eigen::VectorXd& y)
{
	std::cout << t << std::endl;
	return 0;
}

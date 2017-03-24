using namespace lemon;

struct spair{
	int a;
	int b;
	
	spair(int x, int y):a(x),b(y){}
};

bool isOptimal(){
	
}

int lpTask(Matrix &A, std::vector<double> &w, bool logen){
	//Inputs: A matrix and w vector
	
	if (logen) std::cout << "The input matrix: \n" << A << std::endl;
	if (logen) {
		std::cout << "The input w vector: \n" << std::endl;
		for (int i = 0; i < w.size(); i++) 	std::cout << w[i] << "\t";
		std::cout << std::endl;
	}

	std::vector<spair> I,J;
	std::vector<double> v;
	std::vector<std::vector<double> > logm;

	if (logen) std::cout << "Pairs in set I:" << std::endl;
	for (int i = 0; i < 4; i++) {
		v.push_back(log(w[i]));
		std::vector<double> logmrow;
		for (int j = 0; j < 4; j++) {
			logmrow.push_back(log(A.get(i,j)));
			if ((w[i]/w[j] - A.get(i,j)) > pow(10.0,-6.0)) {
				spair tmp(i,j);
				if (logen) std::cout << i << " and " << j << "\n";
				I.push_back(tmp);
			}
		}
		logm.push_back(logmrow);
	}
	if (logen) { 
		std::cout << "The initial v vector (log w): ";
		for (int i = 0; i < v.size(); i++) 	std::cout << v[i] << "\t";
		std::cout << "\n";
	}

	for (int i = 0; i < 4; i++) {
		for (int j = i + 1; j < 4; j++) {
			if (std::abs(w[i]/w[j] - A.get(i,j)) < pow(10.0,-8.0)) {
				spair tmp(i,j);
				J.push_back(tmp);
			}
		}
	}
	if (logen) std::cout << "Sizes of I and J: " << I.size() << " and " << J.size() << "\n";

	Lp LP;
	std::vector<Lp::Col> y; //y=log(x)
	std::vector<Lp::Col> s; //s=log(t)
	if (logen) std::cout << "The LP:\n";
	for (int i = 0; i < 4; i++) { //y-k definialva
		Lp::Col tmp = LP.addCol();
		y.push_back(tmp);
	}
	for (int i = 0; i < I.size(); i++) { //s_ij
		Lp::Col tmp = LP.addCol();
		s.push_back(tmp);
		LP.colLowerBound(tmp, 0);
	}
	for (int i = 0; i < I.size(); i++) {
		Lp::Expr Ex;
		Ex += y[I[i].b];
		Ex -= y[I[i].a];
		LP.addRow(Ex <= (-logm[I[i].a][I[i].b])); //logm a b
		if (logen) std::cout << " y" << I[i].b << " -y" << I[i].a << "<=" << -logm[I[i].a][I[i].b] <<"\n";
		Lp::Expr Ex2;
		Ex2 -= y[I[i].b];
		Ex2 += y[I[i].a];
		Ex2 += s[i];
		LP.addRow(Ex2 <= (v[I[i].a] - v[I[i].b]));
		if (logen) std::cout << "-y" << I[i].b << " +y" << I[i].a << " +s" << i << " <=" << v[I[i].a] <<" - " << v[I[i].b] <<"\n";
	}
	//std::cout << "\n";
	for (int i = 0; i < J.size(); i++) {
		Lp::Expr Ex;
		Ex += y[I[i].a];
		Ex -= y[I[i].b];
		LP.addRow(Ex == logm[I[i].a][I[i].b]);
		if (logen) std::cout << "y" << I[i].a << " -y" << I[i].b << "=" << logm[I[i].a][I[i].b] <<"\n";
	}
	Lp::Expr exy;
	exy += y[0];
	LP.addRow(exy == 0);

	LP.min();
	Lp::Expr e;
	for (int i = 0; i < I.size(); i++) {
		e -= s[i];
	}
	LP.obj(e);
	LP.solve();
	if (logen) std::cout << "The primal type is: " << LP.primalType() << "\n";
	if (LP.primalType() == Lp::OPTIMAL) {
		if (logen){
			std::cout << "Objective function value: " << LP.primal() << std::endl;
			for (int i = 0; i < y.size(); i++) std::cout << "y" << i << ": " << LP.primal(y[i]) << std::endl;
		}
	} else {
		if (logen) {
			std::cout << "Optimal solution not found." << std::endl;
		}
	}
}

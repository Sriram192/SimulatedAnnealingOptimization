using namespace std;
const int maxcolloids = 2000;//what is maxcolloids??
const int maxbin = 200;//what is maxbin?
const int sizeHistVel = 100;
const int bin = 50;
//const double rcut = 2.5;
class Cordinates{
public:
    double x, y;
    void print();
    bool overlap();
}vector;

#define vecsub(a,b,c) c.x = b.x-a.x; c.y= b.y-a.y;
#define dot(a) a.x*a.x + a.y*a.y;
void Cordinates::print(){
	cout << this->x << " " << this->y << " " << &this->x << endl;
}
bool Cordinates::overlap(){
	return ( (this->x > 1000) || (this->y > 1000) );
}
Cordinates colloid[maxcolloids], v[maxcolloids],f[maxcolloids],sumv,sumvn,colloid0[maxcolloids];

class Parameters
{

	double m_a0;
	double m_a1;
	double m_a2;
	double m_a3;
public:
	Parameters(double a0,double a1,double a2, double a3)
	{
	m_a0 = a0;
	m_a1 = a1;
	m_a2 = a2;
	m_a3 = a3;
	}
	
	double get_a0(){ return m_a0;}
	double get_a1(){ return m_a1;}
	double get_a2(){ return m_a2;}
	double get_a3(){ return m_a3;}
		
};

void initialcondition(double &alpha), rdfsample();
double randnum(),RMS(double p),EAV(double q),gauss();
void printout(int t), writeconf(),kinetic(), force(Parameters &input), simulation(Parameters &input), integrate(Parameters &input),rdfsample(), initconf(), anderson(double nu);
double box, alpha, set_tem,rho;
int N, start, ind, baseunit;
double hist[maxbin];//what is hist?
//const int step = 5000;
double U, press, dt=0.005, KE, Temp;
//int index;
double Lindeman();
double t2; //time steps to reach equilibrium
double E_msq ,E_rms,E,Eav,K_msq,Kav,C_v,Sigma_KE,Lin_av,Sigma_E;
double Eq_temp,Eq_press,Eq_pot,Eq_kin;
//double a0,a1,a2,a3;
double tol1;
double MD_NVT(double a0,double a1,double a2,double a3);
double a0_curr,a1_curr,a2_curr,a3_curr;
double a0_best,a1_best,a2_best,a3_best;
double a0_new,a1_new,a2_new,a3_new;
double sim_anneal();
double Ccurr, Cnew, Cbest;
void write_restart_config();
void write_restart_parameters(int iter_index,double &Tcurr, Parameters &Best);
void write_best_parameters(int temp_index,int iter_index,int best_index,double Tcurr,double Cbest,Parameters &Best);
double histVel[sizeHistVel],logHist[sizeHistVel];
double rangeVel,Hfunction;
int min(int x, int y);
double EvalVelDist();
void PrintVelDist();
double EvalVelDist();
int Eq_steps;
double alpha_count(Parameters &Input);
double potential(int i,double x,double y,Parameters &Input);
double totalenergy(Parameters &Input);
double alpha1[bin];
void write_log();
void initialcondition_temp(double &alpha);

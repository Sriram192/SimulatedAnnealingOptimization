//Written by Sriram on 21/11/2018
//MD-NVE code for calculating the melting temperature of the honeycomb potential taken torquato et. al. 2006.
//Also includes a subroutine to calculate the Lindeman Parameter
//Calculation of rdf has to be added according to the discussion with sir on 28/11/2018  //rdf already added.
//Lindeman PArameter has to be printed after every few step and plotted to check the trend and check whether melting is occuring.  //Done successfully
// Start the code with 448 particles to avoid operlap caused by the Initial condition.
//code was corrected by Ramya for overlap of particle images due to PBC.
//Fresh code edited on 04/12/18 for calculation of Cv for calculating the melting temperature.
//start = 0 corresponds to honey comb starting configuration.	
//The formula for C_v has been edited after it was found that the earlier formula was for NVT ensmeble. 
//the formula for C_v given in Frenkel has some mistake refer tildesley for the right formula for NVE
//Refer notes for non dimensionalizing C_v.
//Incorporated calculation of RDF for all the simulations.
//After discussion with sir on 18/12/18, I need to incorporate restart path to get the simulation running from a previously stored configuration(both position and velocity)
//Lot of changes made to sim_anneal() - 23/01/19
//A mistake was found in force calculation and corrected on 28/01/19. this version is the latest to this date(28/01/2019

//-------------------------ImportANT CHANGE-----------------------------------//
//-------------------set temperature changed to 1.0 to accomodate the melting temperature of initial parameters---------//
//--------------------------------------------------------------------------------//
//-----------------22/02/19-------------the cooling schedule changed to a logarithmic schedule apart from changes in the calculation of aceptance ratio-----//
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include<vector>

#include <math.h>
#include <time.h>
#include <ctime>
#include "para.h"
int main (int argc, char *argv[])
{
    srand(time(0));
    if(argc > 1) 
    	{        
	start=atoi(argv[1]);//give 0 for honeycomb
	N = atoi(argv[2]);//give 448
//	alpha = atof(argv[3]);
	}
	
	ind=0;
	set_tem=0.18;

	//box length ;rho is number density;why 1./3. why not 1/3??
                                   /*according to me the function "pow" needs both the arguments
                                   to be of "double" type. therefore we have to enter the aregument
                                   to which the power is raised as a double. Hence 1.0/3.0  *///verify//
//	cout<<"   ind "<<ind<<"   start "<<start<<"   N "<<N<<"   alpha "<<alpha<<"   box "<<box<<"\n";
	
	sim_anneal();


return 0;
}
double MD_NVT(Parameters &input)
{
	alpha = alpha_count(input);
	//function alpha_count calculates the the alpha corresponding to the least potential energy from lattice sums
	cout<<"alpha  "<<alpha<<"\n";
	write_log();
	initialcondition(alpha);
	   
	simulation(input);
	
//	rdfsample();

return Lin_av;
}


void simulation(Parameters &input)
{
//int step = 10000;//my line of code
	cout<<"\nParameters"<<input.get_a0()<<"  "<<input.get_a1()<<"  "<<input.get_a2()<<"  "<<input.get_a3()<<endl;
	E_msq = 0.0,E_rms=0.0,E=0.0,Eav=0.0,K_msq=0.0,Kav=0.0,C_v=0.0,Sigma_KE=0.0,Lin_av=0.0;
	Eq_kin = 0.0, Eq_pot = 0.0, Eq_press = 0.0, Eq_temp = 0.0;
	Eq_steps = 0;
	double tolerance = 1.0;
	int Eq_steps = 50000;
	double E1, E2;
	double H,H1,H2;
//	while(tolerance>0.16)
	for(int i=0;i<Eq_steps;i++)
	{
	
		force(input);
	
		integrate(input);
	
		anderson(3);
	
//		H = EvalVelDist();
//	writeconf();
//		H1 += H;
//		H2 += H*H;
		
		kinetic();
		E = KE+U;
		
		E1 += press;
		E2 += press*press;
		
		
		if(/*Eq_steps*/i%100 ==0) 
		{
	   	 	printout(Eq_steps);
	   	 	cout<<"Eq_time "<</*Eq_steps*/i<<" "<<"\n";
	//    cout<<"particle 12 location "<<&colloid[11].x<<"\n";
		}
		if(/*Eq_steps*/i%100==0)
		{
			E2 /= 100.0;
			E1 /= 100.0;
			tolerance = double(sqrt(E2 - (E1*E1)));//standard deviation
			E2 = 0.0;
			E1 = 0.0;
		}
	}	
//	Eq_steps++;
	
	
	int Pro_steps = 10000;
	
	for(int i=0;i<Pro_steps;i++)
	{
		force(input);
	
		integrate(input);
	
		anderson(3);
	
//	writeconf();
	
		kinetic();	
		
		E = KE+U;
		if(i%100 ==0) 
			{
	 	  	 	printout(Eq_steps);
	 	  	 	cout<<"Pro_steps "<<i<<" "<<"\n";
		
			}
		
		E_msq += E*E;
		Eav += E;
		K_msq += KE*KE;
		Kav += KE;
		Lin_av += Lindeman();
		Eq_kin += KE;
		Eq_pot += U;
		Eq_press += press;
		Eq_temp += Temp;
 		
		
	}
	

	
	Lin_av /= Pro_steps ;
	Kav /= Pro_steps;
	K_msq /= Pro_steps;
	Eav /= Pro_steps;
	E_msq /= Pro_steps;
	Eq_kin /= Pro_steps; //Equilibrium kinetic energy
	Eq_pot /= Pro_steps; //Equilibrium potential energy
	Eq_press /= Pro_steps;
	Eq_temp /= Pro_steps;
	
	
	E_rms = sqrt(E_msq);//my line of code
	Sigma_E = sqrt(E_msq - (Eav*Eav));//standard deviation
	double var;
	var = Kav*Kav;
	
	Sigma_KE = sqrt(K_msq - var);

	double var2;
	var2 = (2.0*(K_msq-var))/(3.0*N*(Temp*Temp));
	
	C_v = double(1.5)/(1-var2);
// 	cout<<"Root mean square kinetic energy "<<K_msq<<"\n";
// 	cout<<"Squate of average KE "<<var<<"\n";
//	cout<<"standard deviation for energy is "<<Sigma_E<<"\n";
//	cout<<"E_RMS value is  "<<E_rms<<"\n";;//my line of code
//	cout<<"specific heat is"<<C_v<<"\n";
	double energy_particle,sigma_energy_particle;
	energy_particle = E_rms/N;
	sigma_energy_particle = Sigma_E/N;
	char STD[100];
	ofstream dev;
    	sprintf( STD, "Equilibrium_values.dat");
    	dev.open (STD, ofstream::out | ofstream::app);
		if (dev.is_open())
    		{
            	dev <<ind<<"  "<< Eq_temp<< "  "<<C_v<<"  "<<Lin_av<<"  "<<Eq_press<<"  "<<Eq_kin<<"  "<<Eq_pot<<"  "<<energy_particle<<" +- "<<sigma_energy_particle<< "\n";
        	}
        	else{ cerr<<"unable to open file for deviation output\n";}
        	ind++;
	dev.close();
return;
}

//----------------------------------------------------------------------
//-----------------------Printing Values--------------------------------
//----------------------------------------------------------------------

void printout(int t)
{
char IntStr[80];
    ofstream ofout;
    sprintf( IntStr, "values%d.dat",ind);
    ofout.open (IntStr, ofstream::out | ofstream::app);
    if (ofout.is_open())
    {
            ofout << ind<<"  "<< t<<"  "<<U<< "  " << KE << "  "<<U+KE<<"  "<<Temp<<"  "<<press<< "\n";
        }
        else {cerr << "unable to open file for config output \n";}
    ofout.close();

//char Lind[100];
//ofstream lin;
//sprintf(Lind, "Lindeman%d.dat",ind);
//lin.open(Lind, ofstream::out| ofstream::app);

//	if(lin.is_open())
//	{
//	lin<< t <<"   "<<Lindeman()<<" \n";
//	}
//	else
//	{
//	cerr<<"unable to open file for Lindeman output\n";
//	}
//lin.close();  

return;
}

//----------------------------------------------------------------------
//--------------------Kinetic Energy&Temp Calculation-------------------
//----------------------------------------------------------------------
void kinetic()
{
	double vsq = 0.0;
	for(int i=0; i<N; i++)
	{
	vsq += v[i].x*v[i].x +v[i].y*v[i].y;//sum of the squares of all the velocities
	}
	KE = vsq / 2.0;
	Temp = vsq/ double(2*N);//This comes from the equipartition theory.

return;
}
//----------------------------------------------------------------------
//-----------------------Initialisation---------------------------------
//----------------------------------------------------------------------
void initialcondition(double &alpha)
{
    if (start ==0){   //honeycomb lattice
	box = pow((double(N)*alpha),1./2.);
	double hexagons = N/2;
	double hexagon_area = double(box*box) / hexagons;

	double baselength = pow(3,0.25) * pow(((2.0*hexagon_area)/9.0),1./2.);
//	cout<<"Base length is  "<<baselength<<"\n";
	int line1 = ceil(box/(cos(0.523) * baselength));
	int line2 = floor(box/(3*baselength));
	 int index = 0;
//	int count = 0;
	
	for(int iy=0;iy<line2;iy++)
		{
		for(int ix=0;ix<=((line1/2)-1);ix++)
			{
			colloid[index].x = ix * 2 * baselength * cos(0.523);
			colloid[index].y = iy * 3 * baselength;
			index+=1;
			}
		}
	for(int iy=0;iy<line2;iy++)
		{
		for(int ix=0;ix<=((line1/2)-1);ix++)
			{
			colloid[index].x = ix * 2 * baselength * cos(0.523);
			colloid[index].y = (iy * 3 * baselength) + baselength;
			index+=1;
			}
		}
	for(int iy=0;iy<line2;iy++)
		{
		for(int ix=0;ix<=((line1/2)-1);ix++)
			{
			colloid[index].x = ix * 2 * baselength * cos(0.523) + (baselength * cos(0.523));
			colloid[index].y = (iy * 3 * baselength) + (1.5*baselength);
			index+=1;
			}
		}
	for(int iy=0;iy<line2;iy++)
		{
		for(int ix=0;ix<=((line1/2)-1);ix++)
			{
			colloid[index].x = ix * 2 * baselength * cos(0.523) + (baselength * cos(0.523));
			colloid[index].y = (iy * 3 * baselength) + (2.5*baselength);
			index+=1;
			}
		}
		
//		cout<<"index"<<index<<"\n";
		N=index;
		box = pow((N*alpha),1./2.);
		for(int i = 0;i<N;i++)
		{
		colloid0[i].x = colloid[i].x;
		colloid0[i].y = colloid[i].y;
		if((colloid0[i].x != colloid[i].x) || (colloid0[i].y != colloid[i].y)) cout<<"doesn't match\n";
		
		}
		//----------------------------
//		Cordinates r; double r2, r1, hbox = box/2.;
//		for(int i=0; i< N-1; i++){
//			for(int j=i+1; j<N; j++){
//				vecsub(colloid[i], colloid[j], r);
//				if (r.x >  hbox) r.x -= box; 
//				else if (r.x < -hbox) r.x += box;
//				if (r.y >  hbox) r.y -= box;
//				else if (r.y < -hbox) r.y += box;
//				r2 = dot(r);
//				r1 = sqrt(r2);
//				if(r1 < 1.0) cout << "overlap " << i << " " << j << "  " << r1 << endl;
//			}
//		}
		
		initconf();	
		}
	
if (start ==1){   //simple cubic
	int index=0;
        int baseunit = ceil(pow(N, 1.0/2.0));/*ceil(x) : Returns the smallest integer that is greater
		                                     than or equal to x (i.e : rounds up the nearest integer).
											 But why use it here?*/
	double lengthcell = box/baseunit;/*the length of the cell or we can say the size of each lattice in the simple cubic structure is
	                                  basically the box length/(number of particles)^(1/3) that is what is done in the previous step 
									  when declaring baseunit = N^(1/3).Basically N^(1/3) ensures that along each length of the box 
									  the particles are arranged in simple cubic fashion*/
	
	  for(int iy=0; iy<baseunit; iy++) {
	     for(int ix=0; ix<baseunit; ix++){
		colloid[index].x=double(ix)*lengthcell;//didnt undertand this?
		colloid[index].y=double(iy)*lengthcell;
		
		index+=1;
             }
	}
	
	for(int i = 0;i<N;i++)
		{
		colloid0[i].x = colloid[i].x;
		colloid0[i].y = colloid[i].y;
		}

initconf();
}






    //---------------------initialize velocities------------------------
    sumv.x = 0.0;
    sumv.y = 0.0;
   
    double vsq = 0.0;
    for(int i=0; i<N; i++){//why for loop initialized from i = 0 rather than i = 1 ?
    v[i].x = (randnum()-0.5);//why randnum() - 0.5??
    v[i].y = (randnum()-0.5);
  
    sumv.x += v[i].x;
    sumv.y += v[i].y;
    
    vsq += v[i].x*v[i].x +v[i].y*v[i].y;   
    }
    cout<<"x momentum"<<" "<<sumv.x<<" "<<"y momentum"<<" "<<sumv.y<<"\n";
    sumv.x /= double(N);
    sumv.y /= double(N);
 
    vsq /= double(N);

    double temp_init = vsq / double(2);
    cout<<"Temperature initialized"<<"  "<<temp_init<<"  "<<"\n";
    double scale = sqrt(2.0*set_tem/vsq);//From where does this scaling come from?
    sumvn.x = 0.0;
    sumvn.y = 0.0;
   
    for(int i=0; i<N; i++){
    	v[i].x = (v[i].x-sumv.x)*scale;
    	v[i].y = (v[i].y-sumv.y)*scale;
    	
        sumvn.x += v[i].x;
        sumvn.y += v[i].y;
        
        vsq += v[i].x*v[i].x +v[i].y*v[i].y;
	}
    vsq /= double(N);
    double temp_init_reset = vsq / double(2.0);
    cout<<"Temperature after reset"<<"  "<<temp_init_reset<<"  "<<"\n";
    cout<<"x momentum"<<" "<<sumvn.x<<" "<<"y momentum"<<" "<<sumvn.y<<"\n";

return;
}
//----------------------------------------------------------------------
double randnum()
{
    return double(rand())/double(RAND_MAX);
}
//----------------------------------------------------------------------
//-----------------------Force Calculation------------------------------
//----------------------------------------------------------------------
void force(Parameters &input)
{
	
double boxinv = 1.0/box;
double boxinv3 = boxinv*boxinv;
double hbox = 0.5*box;
double rcut = 3.0;
double rcutsq = rcut*rcut;
double pi = 3.1415926;
	
for (int i = 0; i<N; i++){//why for loop initialeized from i=0 rather than i=1?
f[i].x=0.0;
f[i].y=0.0;

}
	
int ncut=0;
U = 0;
double W = 0;
Cordinates ri, rij, fij,fi0;
double r1, rij2,sr2,sr3,sr9, sr6, Uij,Wij,fr,UTC,WTC,ideal, r10, r12,r13,r11;
for(int i=0; i<N; i++){
	if(colloid[i].overlap()) {cout << i << endl; colloid[i].print(); cin.ignore();}
	ri.x = colloid[i].x;
	ri.y = colloid[i].y;
	
	fi0.x = f[i].x;
	fi0.y = f[i].y;
	
	for(int j = i+1; j<N; j++){
		rij.x = ri.x - colloid[j].x ;
		rij.y = ri.y - colloid[j].y ;
		
		if (rij.x >  hbox) rij.x -= box; 
		if (rij.x < -hbox) rij.x += box;
		if (rij.y >  hbox) rij.y -= box;
		if (rij.y < -hbox) rij.y += box;
		
		rij2 = rij.x*rij.x + rij.y*rij.y;
		r1 = sqrt(rij2);
		if(r1 < rcut) {
			r12 =1.0/pow(rij2,6);
			r10 = 1.0/pow(rij2,5);
			r13 = 1.0/pow(r1,13);
			r11 = 1.0/pow(r1,11);
				
			U += (5.0*r12) - (input.get_a0()*r10) + input.get_a1()*exp(-input.get_a2()*r1) - 0.4*exp(-40.0*(pow((r1-input.get_a3()),2.0)));
			W += (60.0*r12) - (input.get_a0()*10.0*r10) + (input.get_a1()*input.get_a2()*r1*exp(-input.get_a2()*r1)) - 32*r1*r1*exp(-40.0*(pow((r1-input.get_a3()),2.0))) + 32.0*input.get_a3()*r1*exp(-40.0*(pow((r1-input.get_a3()),2.0)));
			
			fr = (60.0*r13) - (input.get_a0()*10.0*r11) + (input.get_a1()*input.get_a2()*exp(-input.get_a2()*r1)) - 32.0*r1*exp(-40.0*(pow((r1-input.get_a3()),2.0))) + 32.0*input.get_a3()*exp(-40.0*(pow((r1-input.get_a3()),2.0)));
			fij.x = (1.0/r1)*fr*rij.x;
			fij.y = (1.0/r1)*fr*rij.y;
			
			fi0.x += fij.x;
			fi0.y += fij.y;
				
			f[j].x -= fij.x; // force on particle j is equal and opposite to the force on the i particle for that particular ij pair
			f[j].y -= fij.y;
				
			ncut += 1;
		}
	}
	f[i].x = fi0.x;
	f[i].y = fi0.y;
	
}
//sr2 = 1.0/rcutsq;
//sr6 = sr2*sr2*sr2;
//sr3 = sr2/rcut;
//sr9 = sr3*sr3*sr3;
//Uij = sr6*(sr6 - 1.0);
//U -= double(ncut)*Uij;

//for(int i=0; i<N; i++){
//f[i].x *= 24.0;
//f[i].y *= 24.0;

//}

//W = W/3.0;
//U *= 4.0;
//UTC = 8./9.*pi*(1/alpha)*double(N)*(sr9 - 3.0*sr3);
//WTC = 16./9.*pi*(1/alpha)*(1/alpha)*double(N)*(2.*sr9 - 3.*sr3);
//U += UTC;
//W += WTC;
ideal = double(N)*boxinv3*set_tem;
press = ideal+W*boxinv3;
//cout<<"UTC     "<<UTC<<" "<<"\n";

return;
}
//----------------------------------------------------------------------
//---------------------------Integration--------------------------------
//
//----------------------------------------------------------------------
void integrate(Parameters &input)
{
double dt2=dt/2.0;
double dtsq2 = dt*dt2;
double hbox = 0.5*box;


for(int i=0; i<N; i++) {
	colloid[i].x += (dt*v[i].x + dtsq2*f[i].x);
	colloid[i].y += (dt*v[i].y + dtsq2*f[i].y);
		if(colloid[i].overlap()) {cout << i << endl; colloid[i].print(); v[i].print(); f[i].print(); cin.ignore();}
	
	if(colloid[i].x >= box) colloid[i].x -= box;
	if(colloid[i].x <= 0.0) colloid[i].x += box;
 	if(colloid[i].y >= box) colloid[i].y -= box;
        if(colloid[i].y <= 0.0) colloid[i].y += box;
 	

	v[i].x += (dt2* f[i].x);
	v[i].y += (dt2* f[i].y);
	
}

force(input);

for(int i=0; i<N; i++) {
	v[i].x += (dt2* f[i].x);
	v[i].y += (dt2* f[i].y);
	
}

return;
}

//----------------------------------------------------------------------
//---------------------Printing coordinates-----------------------------
//----------------------------------------------------------------------
//void writeconf()
//{
//    char IntStr[80];
//    ofstream of;
//    sprintf( IntStr, "config%d.xyz",ind);
//    of.open (IntStr, ofstream::out | ofstream::app);
//    if (of.is_open())
//    {
//        of << N<< "\n";
//        of << "\n";
//	if ((start==0) || (start ==1)){
//        for (int i=0; i<N; i++)
//        {
//            of << 'h'<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< 0.0 << "\n";
//        }    
//    } 
//	else if (start ==3) {
//	 for (int i=0; i<N; i++)
//        {
//	   if (i%4 ==0) {
//            of <<"he"<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< 0.0 << "\n";
//	   }else {of << 'h'<<"  "<<colloid[i].x << "  " << colloid[i].y << "  "<< 0.0 << "\n";}
//	}
//        }

//        }
//	else {cerr << "unable to open file for config output \n";}
//    of.close();
//    
//    return;
//}
//---------------------------------------------------------------------------------------//
//my lines of code---------------------
//double RMS(double p)
//{
//	double E1;
//	int step =15000;
//	E1 = (p*p)/step;
//	
//	
//	return(E1);
//}
//double EAV(double q)
//{
//	double E2;
//	int step = 15000;
//	E2 = q/step;
//	
//	return(E2);
//}
void initconf()
{
char conf[100];
ofstream val;
sprintf(conf,"Initialconf.xyz");
val.open(conf,ofstream::out | ofstream::trunc);

	if(val.is_open())
	{
	val<<N+4<<"\n";
	val<<"\n";
	
	val<<'c'<<"  "<<0.0<<"  "<<0.0<<"  "<<0.0<<"\n"; 
	val<<'c'<<"  "<<0.0<<"  "<<box<<"  "<<0.0<<"\n"; 
	val<<'c'<<"  "<<box<<"  "<<0.0<<"  "<<0.0<<"\n"; 
	val<<'c'<<"  "<<box<<"  "<<box<<"  "<<0.0<<"\n";  
	for(int i=0; i<N; i++)
		{
		val<<'h'<<"  "<<colloid[i].x<<"  "<<colloid[i].y<<"  "<<0.0<<"\n";
		}
	}
	else
	{
	cerr<<"unable to open outpur file \n";
	}
	val.close();
	
}
//-----------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------Lindeman  Parameter-----------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------------------//

double Lindeman()
{
//calculating the lindmann parameter for the MD runs
	double r,val1, val2,Lin,rij2;
	Cordinates ri;
	val1 = 0;
	val2 = 0;
	for(int i=0; i<N; i++)
	{
	vecsub(colloid0[i],colloid[i],ri)
//	ri.x = colloid[i].x - colloid0[i].x;
//	ri.y = colloid[i].y - colloid0[i].y;
//	
	rij2 = ri.x*ri.x + ri.y*ri.y;
	r = sqrt(rij2);
	
	val1 += rij2;
	val2 += r/N;
	}
	val1 /= N;
	val2 *= val2;
	
	Lin = sqrt(val1-val2);
	
//	cout<<"Lindeman Parameter  "<<Lin<<"\n";

return Lin;	
}
//----------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------Radial distribution Function---------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------// 
//void rdfsample()
//{
//	double delr = box/double(2*maxbin);//why 2*maxbin??
//	double pi = 3.1415926;
//	double boxinv = 1.0/box;//why are we calaculating this value?
//	double hbox = 0.5*box;//half box length
//	double rij2, r1;//rij2 is rij^square and r1 is the magnitude of the distance between the particles
//	int nsample = 0, bin;
//	Cordinates rij;
//	char IntStr[80];
//	ofstream ofrdf;
//	sprintf( IntStr, "rdf%d.dat",ind);
//    	ofrdf.open (IntStr, ofstream::out | ofstream::app);
//	for(int i=0; i<maxbin; i++)
//	{
//	hist[i]=0;//what is hist?
//	}

//	for (int t=0; t<5000; t++)
//	{
//	force();
//	integrate();
//	anderson(3);
//		if(t%1 ==0) 
//		{
//		cout<<t<<"\n";
//		nsample += 1;
//			for (int i = 0; i<N; i++) 
//				{
//    				for(int j=i+1; j<N; j++) 
//    					{
//       					rij.x = colloid[i].x - colloid[j].x ;
//       					rij.y = colloid[i].y - colloid[j].y ;
//     
//       
//					if (rij.x >  hbox) rij.x -= box;//Periodic boundary condition
//       					if (rij.x < -hbox) rij.x += box;//Periodic boundary condition
//				       	if (rij.y >  hbox) rij.y -= box;//Periodic boundary condition
//       					if (rij.y < -hbox) rij.y += box;//Periodic boundary condition
//       					
//       					rij2 = rij.x*rij.x + rij.y*rij.y;
//       					r1 = sqrt(rij2);
//       /*Instead of calculating square root of "r" outside the if condition
//        *we can check whether r1^2 < hbox^2 and then if the condition satisfies
//        *we calcuate the square root of "r". this will be computationally 
//		*less intensive*/
//       						if(r1 <= hbox)
//       						{
//        					bin = int(r1/delr);//Address of the particle i.e the bin number to which the particle j belongs wrt the ith partice. 
//						hist[bin] += 2.0;//this is because there are pairs. jthe particle belongs to bin and vice versa fo jthe particle
//						} 
//					}
//				}
//		}
//	}
//	rho = 1.0/alpha;
//	double idcon = 2.0*pi*rho;
//	double rlow, rup, ideal, rbin;
//		if(ofrdf.is_open()) 
//		{
//			for(int i=0; i<maxbin; i++)
//			{
//			rlow = double(i)*delr;
//			rup = rlow + delr;
//			ideal = idcon*(rup*rup - rlow*rlow);
//			hist[i] = hist[i]/double(N*nsample)/ ideal;
//			rbin = rlow + delr/2.0;
//            	ofrdf << rbin <<"  "<<hist[i]<< "\n";
//			}
//		}
//		else {cerr << "unable to open file for rdf output \n";}
//    ofrdf.close();
//return;
//}
//----------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------Simulated Annealing optimization--------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------//
double sim_anneal()
{
//	double a0_curr,a1_curr,a2_curr,a3_curr;
//	double a0_best,a1_best,a2_best,a3_best;
//	double a0_new,a1_new,a2_new,a3_new;
//	Parameters Curr;
//	Parameters New;
//	Parameters Best;
	
	//Initiliazing randum values for potential parameters.
	a0_curr = 6.5;
	a1_curr = 18.5;
	a2_curr = 2.45;
	a3_curr = 1.83;
	Parameters Curr(a0_curr,a1_curr,a2_curr,a3_curr);
	
	//Getting best parameter values
	a0_best = Curr.get_a0();
	a1_best = Curr.get_a1();
	a2_best = Curr.get_a2();
	a3_best = Curr.get_a3();
	Parameters Best(a0_best,a1_best,a2_best,a3_best);
	
	//Initial temperature for the simulated annealing optimization loop
	double Tcurr;
	Tcurr = 7.0; //arbitrary choice for the temperature. The maximum of the cost function gives us the temperature generally. Not used here.
	
	
	int temp_index = 0;
	int best_index = 1;
	bool update;
	
	double Ccurr=0.0, Cnew = 0.0;
//	std::vector<double> Cbest;
//	double Cbest;
//	cout<<"size of Cbest"<<Cbest.size()<<"\n";
//	tol1 = 5.0;
	int iter_index = 0;
	int acceptance_index = 0;
	double acc_ratio;
	while(Tcurr>2.0) 
	{	
	tol1 = 5.0;
		while(tol1 >= 5.0)	
		{
			Cnew = 0.0; Ccurr = 0.0; Cbest = 0.0;
			//assingning new values to parameters
			a0_new = Curr.get_a0() + 0.1*(randnum()-0.5);
			a1_new = Curr.get_a1() + 0.1*(randnum()-0.5);
			a2_new = Curr.get_a2() + 0.1*(randnum()-0.5);
			a3_new = Curr.get_a3() + 0.1*(randnum()-0.5);
			Parameters New(a0_new,a1_new,a2_new,a3_new);
			
				for(int i=0;i<3;i++)
				{
				Ccurr += MD_NVT(Curr);
				}
			cout<<"\n------------------------------------Curr over-------------------------------\n";
				for(int i=0;i<3;i++)
				{
				Cnew += MD_NVT(New);//new parameters
				}
			cout<<"\n------------------------------------New over-------------------------------\n";
			Ccurr /= 3.0;
			Cnew /= 3.0;
				cout<<"Cnew	"<<Cnew<<"	"<<"Ccurr	"<<Ccurr<<"\n";

				if(Cnew<Ccurr)
				{
				cout<<"\n Inner if loop \n";	
				tol1 = ((Ccurr-Cnew)/Ccurr) * 100.0;
				//updating the current parameters
				Curr = New;
//				
					for(int i=0;i<3;i++)
					{                  
					Cbest += MD_NVT(Best);
					}
					Cbest /= 3.0;
					cout<<"Cbest	"<<Cbest<<"	"<<"Cnew	"<<Cnew<<"\n";
				cout<<"\n------------------------------------Best over-------------------------------\n";
						if(Cnew<Cbest)
						{
						cout<<"\n Inner of inner loop\n";
						Best = New;
						best_index++;
//						a0_best = New.get_a0();
//						a1_best = New.get_a1();
//						a2_best = New.get_a2();
//						a3_best = New.get_a3();
//						Parameters Best(a0_best,a1_best,a2_best,a3_best);
						write_best_parameters(temp_index,iter_index,best_index,Tcurr,Cbest,Best);
						
//						update = true;					
						}
//						else
//						{
//						update = false;
//						int x;
//						cout<<"enter x \n";
//						cin>>x;
//						cout<<"\n";
//						}
				
				
//				acceptance_index++;
				}
				else if(randnum()<(exp(-(Cnew-Ccurr)/Tcurr)))  //simulated annealing acceptance crieteria
					{
					cout<<"\n Inner else if loop\n";
					tol1 = ((Cnew - Ccurr)/Ccurr)*100;
					Curr = New;
					acceptance_index++;
					}
				
				iter_index++;
				
		cout<<"\niter_index "<<iter_index<<"\n"; 
		
//		if(update==true)
//			{
//				
//			}
			
			char OutStr1[100];
			ofstream out2;
			sprintf(OutStr1,"Potential_parameters.dat"); 
			out2.open(OutStr1,ofstream::out | ofstream::app);
			
				if(!out2) 
				{
				cerr<<"Unable to open file to write potential parameters\n";
				}
			
				if(out2.is_open())
				{
					out2<<iter_index-1<<"  "<<"  "<<Tcurr<<"  "<<Cnew<<"  "<<New.get_a0()<<"  "<<New.get_a1()<<"  "<<New.get_a2()<<"  "<<New.get_a3()<<"\n";
				}
				out2.close();
		
		
		
//		iter_index++;		
		cout<<"\n inner loop \n";
		acc_ratio = double(acceptance_index/iter_index);
		cout<<"acceptance ratio--------------> "<<acc_ratio<<"\n";
		write_restart_config();
		write_restart_parameters(iter_index,Tcurr,Best);
		}
		
		
		Tcurr /= log(3.0);
		temp_index++;
		cout<<"\n outer loop\n";
	}
	
return 0; //make sure that after the optimization Xbest is returned or printed out.
}

double gauss()//what does this function do?
{
    double pr=2,v1,v2,l1;
    do{
        v1 = 2.*randnum()-1.0;
        v2 = 2.*randnum()-1.0;
        pr=v1*v1+v2*v2;
    }while(pr >=1.0);
        l1 = v1*sqrt(-2*log(pr)/(pr));
    return sqrt(set_tem)*l1;
}


void anderson(double nu)
{
    int i;
    double T_req=set_tem,ran[N];
    double u1_ran, u2_ran, u3_ran, a1, a2, a3;
    float n1_ran, n2_ran, n3_ran, pi=3.1415;
    for(i=0;i<N;i++)
    {
        if(randnum()<=nu*dt)
        {
            v[i].x=gauss();
            v[i].y=gauss();
//            v[i].z=gauss();
            //cout<<v[i].x<<"\t"<<v[i].y<<"\t"<<v[i].z<<"\n";

        }                                                      
    }
}	
void write_restart_config()
{
char rest[100];
ofstream op;
sprintf(rest,"config_restart.xyz");
op.open(rest,ofstream::out | ofstream::trunc);
	if(!(op.is_open())){cerr<<"Unable to open file for config restart output";}
	op<< N<<"\n";
	op<<"\n";
	for(int i=0;i<N;i++)
	{
		op<<colloid[i].x<<"  "<<colloid[i].y<<"\n";
	}
op.close();
}				
void write_restart_parameters(int iter_index,double &Tcurr, Parameters &Best)
{
char rest[100];
ofstream op;
sprintf(rest,"parameter_restart.dat");
op.open(rest,ofstream::out | ofstream::trunc);
	if(!(op.is_open())){cerr<<"Unable to open file for parameter restart output";}
//	op<< N<<"\n";
//	op<<"\n";
		op<<Tcurr<<"  "<<Best.get_a0()<<"  "<<Best.get_a1()<<"  "<<Best.get_a2()<<"  "<<Best.get_a3()<<"\n";  
op.close();
}
//----------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------H - Function-------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------//
double EvalVelDist()
{
	double deltaV,histSum;
	int j, n;
	double rangeVel = 5;
	double vdot;
	
	double Hfunction;
	
	
	for(int j=0;j<sizeHistVel;j++){
	
	histVel[j] = 0.0;
}
   
	deltaV = rangeVel/sizeHistVel;
	
	for(int i=0;i<N;i++){
	
	vdot = (v[i].x*v[i].x) + (v[i].y*v[i].y);
	vdot = sqrt(vdot);
	int k = vdot/deltaV;
	
	int c = min(k,sizeHistVel-1);
	
	histVel[c] += 1; 
	
	
	}

		histSum = 0.0;
		
			
	
		
		for(int j=0;j<sizeHistVel;j++){
			if(histVel[j]>0){
			
			logHist[j] = log (histVel[j]);
			
			histSum += (logHist[j] * histVel[j] * deltaV);
			}
			
		}
	
//		Hfunction = 0.5*histSum;
//		//cout<<Hfunction;
//		char IntStr[80];
//    ofstream ofout;
//    sprintf( IntStr, "Hfun%d.dat",ind);
//    ofout.open (IntStr, ofstream::out | ofstream::app);
//    
//            ofout << Eq_steps<<"  "<<Hfunction<< "\n";
//    ofout.close()    ;
		
//			PrintVelDist();
		

	
return Hfunction;
}

int min(int x, int y)
{
	int a;
	if(x>y){
	
	a=y;
}
	else{
		a=x;
	}
	return a;
}
//void PrintVelDist()
//{
//	char IntStr[80];
//	double vBin;
//	ofstream ofvel;
//    sprintf( IntStr, "veld%d.dat",ind);
//    ofvel.open (IntStr, ofstream::out | ofstream::app);
//    for(int n=0;n<sizeHistVel;n++){
//    	//vBin = (n+0.5) * (rangeVel / sizeHistVel);
//    	//cout<<vBin;
//    	ofvel << n <<"  "<<histVel[n]<< "\n";
//	ofvel.close();
//	}
//	return;
//}
double alpha_count(Parameters &Input)
{

	char pot_energy[100];
	ofstream pot;
    	sprintf( pot_energy, "energy_per_particle.dat");
    	pot.open (pot_energy, ofstream::out | ofstream::app);
	int a_count = 0;
	alpha1[0] = 0.7;
	double least_energy;
	double least = 100.0;
	int least_index = 0;
		for(int i=1; i<50;i++)
		{
		alpha1[i] = alpha1[i-1] + 0.05;
		a_count += 1;
		}
	
	
		for(int i=0;i<a_count;i++)
		{
		initialcondition_temp(alpha1[i]);
		//cout<<"after initial \n";
		least_energy/*[i]*/ = double(totalenergy(Input))/N;
			if (pot.is_open())
			{
			pot<< alpha1[i] <<"	"<<double(totalenergy(Input))/N<<"\n";
			}
			else{cerr<<"Unable to open file for writing output";}
		if(least_energy/*[i]*/<least)
			{
			least = least_energy/*[i]*/;
		
			least_index = i;
			}
		}
	pot.close();	
	
	

return alpha1[least_index];
}
double totalenergy(Parameters &Input)
{
	double pot =0;
		for(int i=0; i<N;i++)
		{
		pot += potential(i,colloid[i].x,colloid[i].y,Input);
		}
//double r3 = 1.0/(cutoff*cutoff*cutoff);
//double Utail = 8*3.14*density*r3*number_of_colloids*(((1.0/9.0)*(pow(r3,2))) - (1.0/3.0));
return (double(pot/2));// + Utail);
}
double potential(int i,double x,double y,Parameters &Input)
{
	double xij, yij;
//In actual simulation, get the cutoff radius as an input.
	double cutoff = 2.5, cutoff2;
	double r2,r,U=0.0,r12,r10;

		for(int j=0; j<N; j++)
		{
			if(j != i)
			{
			xij = double(x-colloid[j].x);
			yij = double(y-colloid[j].y);
//		zij = double(z-colloid_z[j]);
//minimum image convention	
//		xij = periodic(xij);
//		yij = periodic(yij);
//		zij = periodic(zij);

//in actual simulation, check for periodic BC
			r2 = double((xij*xij) + (yij*yij));
		
			cutoff2 = double(cutoff*cutoff);
	
			if(r2 < cutoff2)
			{
			r = sqrt(r2);//defining radial distance from the particle
			r12 =1.0/pow(r2,6);
			r10 = 1.0/pow(r2,5);
			U += (5.0*r12) - (Input.get_a0()*r10) + Input.get_a1()*exp(-Input.get_a2()*r) - 0.4*exp(-40.0*(pow((r-Input.get_a3()),2.0)));
			}
			}
		}

return U;
}			
void write_log()
{
char logfile[100];
	ofstream log;
	sprintf(logfile,"log.dat");
	log.open(logfile,ofstream::out | ofstream::trunc);
	if(log.is_open())
	{
	log<<"Number of particles  "<<N<<"\n";
	log<<"Area fraction(alpha)  "<<alpha<<"\n";
	if(start==0){log<<"Initial config (start)  "<<start<<"  Honeycomb"<<"\n";}
	if(start==1){log<<"Initial config (start)  "<<start<<"  Simple cubic"<<"\n";}
	log<<"Initial temperature  "<< set_tem<<"\n";
	}else{cerr<<"unable to open log file\n";}
	log.close();
}	
void initialcondition_temp(double &alpha)
{
       //honeycomb lattice
	box = sqrt(N*alpha);
	double hexagons = N/2;
	double hexagon_area = double(box*box) / hexagons;

	double baselength = pow(3,0.25) * pow(((2.0*hexagon_area)/9.0),1./2.);
	
	int line1 = ceil(box/(cos(0.523) * baselength));
	int line2 = floor(box/(3*baselength));
	int index = 0;
//	int count = 0;
	
	for(int iy=0;iy<line2;iy++)
		{
		for(int ix=0;ix<=((line1/2)-1);ix++)
			{
			colloid[index].x = ix * 2 * baselength * cos(0.523);
			colloid[index].y = iy * 3 * baselength;
			index+=1;
			}
		}
	for(int iy=0;iy<line2;iy++)
		{
		for(int ix=0;ix<=((line1/2)-1);ix++)
			{
			colloid[index].x = ix * 2 * baselength * cos(0.523);
			colloid[index].y = (iy * 3 * baselength) + baselength;
			index+=1;
			}
		}
	for(int iy=0;iy<line2;iy++)
		{
		for(int ix=0;ix<=((line1/2)-1);ix++)
			{
			colloid[index].x = ix * 2 * baselength * cos(0.523) + (baselength * cos(0.523));
			colloid[index].y = (iy * 3 * baselength) + (1.5*baselength);
			index+=1;
			}
		}
	for(int iy=0;iy<line2;iy++)
		{
		for(int ix=0;ix<=((line1/2)-1);ix++)
			{
			colloid[index].x = ix * 2 * baselength * cos(0.523) + (baselength * cos(0.523));
			colloid[index].y = (iy * 3 * baselength) + (2.5*baselength);
			index+=1;
			}
		}
		
		N=index;
		box = pow((N*alpha),1./2.);
		
		//----------------------------
//		Cordinates r; double r2, r1, hbox = box/2.;
//		for(int i=0; i< N-1; i++)
//		{
//			for(int j=i+1; j<N; j++)
//			{
//				vecsub(colloid[i], colloid[j], r);
//				if (r.x >  hbox) r.x -= box; 
//				else if (r.x < -hbox) r.x += box;
//				if (r.y >  hbox) r.y -= box;
//				else if (r.y < -hbox) r.y += box;
//				r2 = dot(r);
//				r1 = sqrt(r2);
//				if(r1 < 1.0) cout << "overlap " << i << " " << j << "  " << r1 << endl;
//			}
//		}
		
		
}
void write_best_parameters(int temp_index,int iter_index,int best_index,double Tcurr,double Cbest,Parameters &Best)
{
cout<<"Printing values\n";
			char OutStr[100];
			ofstream out1;
			sprintf(OutStr,"Best_Parameters%d.dat",temp_index); 
			out1.open(OutStr,ofstream::out | ofstream::app);
			
				if(!out1) 
				{
				cerr<<"Unable to open file to write potential parameters\n";
				}
			
				if(out1.is_open())
				{
					out1<<iter_index<<"  "<<best_index<<"  "<<Tcurr<<"  "<<Cbest<<"  "<<Best.get_a0()<<"  "<<Best.get_a1()<<"  "<<Best.get_a2()<<"  "<<Best.get_a3()<<"\n";
				}
				out1.close();
			cout<<"Exiting Printing\n";
}

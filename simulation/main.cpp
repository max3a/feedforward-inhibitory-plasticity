#include "run.h"

int main(int argc,char*argv[])
{
	double theta;
	double gamma;
	int pulse_number;	
	int count;
	int t_record;
	int LTP;
	
	
	////////////////setting
	t_record=1;//record synaptic current and membrane potential every t_record ms
	LTP=1; //feedforward plasticity :0 off, 1 on
	
	count=26380; //numbering of simulation
	
	//parameter of rescue
	theta=8; //theta frequency
	gamma=40; //gamma frequency
	pulse_number=1; //number of pulse applied within one theta period, if LTP=1, pulse_number will be set to 0	
	///////////////		
	

	evolution(count,time(NULL),theta,gamma,pulse_number,t_record,LTP);
}

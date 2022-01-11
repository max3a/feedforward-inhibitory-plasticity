# feedforward-inhibitory-plasticity
C++ code for the paper "Rational designing of oscillatory rhythmicity for memory rescue in plasticity-impaired learning networks". You need a gnu complier, and RandNum.h and run.h provided whitin this repository to run the program.

In main.cpp, you can decide 
    -the recording time period of synaptic current, t_record. 
    -Able/disable feedforward E->I plasticity, LTP.
    -the theta, and gamma frequency and number of pulses applied within one theta period during rescue when feedforward E->I plasticity is not disabled, theta, gamma,                  pulse_number. 
    -the numbering of the simulation, count

To run the simulation, you can run cmd in window, go to the directory folder and type mingw32-make. It will generate an exe, named "main.exe".

By running "main.exe", 9 txt files will be generated, which are
    -E=3.00I=8.00P=0.20_spike_{count}\_0_0.00.txt

sprintf(filename_2,"E=%.2fI=%.2fP=%.2f_spike_%d_0_%.2f.txt",x,y,link_density,count,yy);
	sprintf(filename_3,"test_%d_RE_AMPA.txt",count);
	sprintf(filename_3,"test_%d_RE_NMDA.txt",count);
	sprintf(filename_3,"test_%d_RE_GABA.txt",count);
	sprintf(filename_3,"test_%d_FF_AMPA.txt",count);
	sprintf(filename_3,"test_%d_FF_NMDA.txt",count);
	sprintf(filename_3,"test_%d_FF_GABA.txt",count);
	sprintf(filename_3,"test_%d_V.txt",count);
    sprintf(filename_5,"E=%.2fI=%.2fP=%.2f_W_%d_%.0f_%.2f.txt",x,y,link_density,count,t*delta_t,yy);

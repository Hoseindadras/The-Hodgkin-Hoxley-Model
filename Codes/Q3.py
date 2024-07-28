from brian2 import *
import matplotlib.pyplot as plt

cohs = [20,17,15,12.5,5,2,1,0,-1,-2,-5,-12,-17,-20]
for cohid,coh in enumerate(cohs):
   coh = -12.8 
   sigma = 4.0 * Hz 
   mu0 = 40.0 * Hz  
   mu1 = 40.0 * Hz  
   stim_interval = 50.0 * ms  
   stim_on = 1000 * ms  
   stim_off = 3000 * ms  
   runtime = 4000 * ms  

   N_ext = 1000 
   rate_ext_E = 2400 * Hz / N_ext  
   rate_ext_I = 2400 * Hz / N_ext  

   N = 2000  
   f_inh = 0.2  
   NE = int(N * (1.0 - f_inh))  
   NI = int(N * f_inh)  
   fE = 0.15 
   subN = int(fE * NE)  

   El = -70.0 * mV  
   Vt = -50.0 * mV  
   Vr = -55.0 * mV 
   CmE = 0.5 * nF  
   CmI = 0.2 * nF 
   gLeakE = 25.0 * nS 
   gLeakI = 20.0 * nS  
   refE = 2.0 * ms  
   refI = 1.0 * ms  

   V_E = 0. * mV 
   V_I = -70. * mV  
   tau_AMPA = 2.0 * ms 
   tau_NMDA_rise = 2.0 * ms  
   tau_NMDA_decay = 100.0 * ms  
   tau_GABA = 5.0 * ms  
   alpha = 0.5 * kHz 
   C = 1 * mmole  

   gextE = 2.1 * nS  
   gextI = 1.62 * nS 
   gEEA = 0.05 * nS / NE * 1600 
   gEIA = 0.04 * nS / NE * 1600 
   gEEN = 0.165 * nS / NE * 1600 
   gEIN = 0.13 * nS / NE * 1600  
   gIE = 1.3 * nS / NI * 400 
   gII = 1.0 * nS / NI * 400  


   Jp = 1.7  
   Jm = 1.0 - fE * (Jp - 1.0) / (1.0 - fE)

  
   eqsE = """
      label : integer (constant)  # label for decision encoding populations
      dV/dt = (- gLeakE * (V - El) - I_AMPA - I_NMDA - I_GABA - I_AMPA_ext + I_input) / CmE : volt (unless refractory)

      I_AMPA = s_AMPA * (V - V_E) : amp
      ds_AMPA / dt = - s_AMPA / tau_AMPA : siemens

      I_NMDA = gEEN * s_NMDA_tot * (V - V_E) / ( 1 + exp(-0.062 * V/mvolt) * (C/mmole / 3.57) ) : amp
      s_NMDA_tot : 1

      I_GABA = s_GABA * (V - V_I) : amp
      ds_GABA / dt = - s_GABA / tau_GABA : siemens

      I_AMPA_ext = s_AMPA_ext * (V - V_E) : amp
      ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : siemens

      I_input : amp

      ds_NMDA / dt = - s_NMDA / tau_NMDA_decay + alpha * x * (1 - s_NMDA) : 1
      dx / dt = - x / tau_NMDA_rise : 1
   """

   eqsI = """
      dV/dt = (- gLeakI * (V - El) - I_AMPA - I_NMDA - I_GABA - I_AMPA_ext) / CmI : volt (unless refractory)

      I_AMPA = s_AMPA * (V - V_E) : amp
      ds_AMPA / dt = - s_AMPA / tau_AMPA : siemens

      I_NMDA = gEIN * s_NMDA_tot * (V - V_E) / ( 1 + exp(-0.062 * V/mvolt) * (C/mmole / 3.57) ): amp
      s_NMDA_tot : 1

      I_GABA = s_GABA * (V - V_I) : amp
      ds_GABA / dt = - s_GABA / tau_GABA : siemens

      I_AMPA_ext = s_AMPA_ext * (V - V_E) : amp
      ds_AMPA_ext / dt = - s_AMPA_ext / tau_AMPA : siemens
   """

   popE = NeuronGroup(NE, model=eqsE, threshold='V > Vt', reset='V = Vr', refractory=refE, method='euler', name='popE')
   popI = NeuronGroup(NI, model=eqsI, threshold='V > Vt', reset='V = Vr', refractory=refI, method='euler', name='popI')
   popE1 = popE[:subN]
   popE2 = popE[subN:2 * subN]
   popE3 = popE[2 * subN:]
   popE1.label = 0
   popE2.label = 1
   popE3.label = 2

   C_EE_AMPA = Synapses(popE, popE, 'w : siemens', on_pre='s_AMPA += w', delay=0.5 * ms, method='euler', name='C_EE_AMPA')
   C_EE_AMPA.connect()
   C_EE_AMPA.w[:] = gEEA
   C_EE_AMPA.w["label_pre == label_post and label_pre < 2"] = gEEA*Jp
   C_EE_AMPA.w["label_pre != label_post and label_post < 2"] = gEEA*Jm
 
   C_EI_AMPA = Synapses(popE, popI, on_pre='s_AMPA += gEIA', delay=0.5 * ms, method='euler', name='C_EI_AMPA')
   C_EI_AMPA.connect()

   C_EE_NMDA = Synapses(popE, popE, on_pre='x_pre += 1', delay=0.5 * ms, method='euler', name='C_EE_NMDA')
   C_EE_NMDA.connect(j='i')

   NMDA_sum_group = NeuronGroup(3, 's : 1', name='NMDA_sum_group')

   NMDA_sum = Synapses(popE, NMDA_sum_group, 's_post = s_NMDA_pre : 1 (summed)', name='NMDA_sum')
   NMDA_sum.connect(j='label_pre')

   NMDA_set_total_E = Synapses(NMDA_sum_group, popE,
                              '''w : 1 (constant)
                                 s_NMDA_tot_post = w*s_pre : 1 (summed)''', name='NMDA_set_total_E')
   NMDA_set_total_E.connect()
   NMDA_set_total_E.w = 1
   NMDA_set_total_E.w["i == label_post and label_post < 2"] = Jp
   NMDA_set_total_E.w["i != label_post and label_post < 2"] = Jm

   NMDA_set_total_I = Synapses(NMDA_sum_group, popI,
                              '''s_NMDA_tot_post = s_pre : 1 (summed)''', name='NMDA_set_total_I')
   NMDA_set_total_I.connect()

   C_IE = Synapses(popI, popE, on_pre='s_GABA += gIE', delay=0.5 * ms, method='euler', name='C_IE')
   C_IE.connect()

   C_II = Synapses(popI, popI, on_pre='s_GABA += gII', delay=0.5 * ms, method='euler', name='C_II')
   C_II.connect()

   extinputE = PoissonInput(popE, 's_AMPA_ext', N_ext, rate_ext_E, gextE)
   extinputI = PoissonInput(popI, 's_AMPA_ext', N_ext, rate_ext_I, gextI)

   stiminputE1 = PoissonGroup(subN, rates=0*Hz, name='stiminputE1')
   stiminputE2 = PoissonGroup(subN, rates=0*Hz, name='stiminputE2')
   stiminputE1.run_regularly("rates = int(t > stim_on and t < stim_off) * (mu0 + coh / 100.0 * mu1 + sigma*randn())", dt=stim_interval)
   stiminputE2.run_regularly("rates = int(t > stim_on and t < stim_off) * (mu0 - coh / 100.0 * mu1 + sigma*randn())", dt=stim_interval)
   C_stimE1 = Synapses(stiminputE1, popE1, on_pre='s_AMPA_ext += gextE', name='C_stimE1')
   C_stimE1.connect(j='i')
   C_stimE2 = Synapses(stiminputE2, popE2, on_pre='s_AMPA_ext += gextE', name='C_stimE2')
   C_stimE2.connect(j='i')


   popE.s_NMDA_tot = tau_NMDA_decay * 10 * Hz * 0.2
   popI.s_NMDA_tot = tau_NMDA_decay * 10 * Hz * 0.2
   popE.V = Vt - 2 * mV
   popI.V = Vt - 2 * mV

   SME1 = SpikeMonitor(popE1, record=True)
   SME2 = SpikeMonitor(popE2, record=True)

   R1 = PopulationRateMonitor(popE1)
   R2 = PopulationRateMonitor(popE2)

   E1 = StateMonitor(stiminputE1, 'rates', record=0, dt=1*ms)
   E2 = StateMonitor(stiminputE2, 'rates', record=0, dt=1*ms)

   run(runtime, report='stdout', profile=True)
   print(profiling_summary())

   fig, axs = plt.subplots(4, 1, sharex=True, layout='constrained', gridspec_kw={'height_ratios': [2, 2, 2, 1]})
   axs[0].plot(SME1.t / ms, SME1.i, '.', markersize=2, color='darkred')
   axs[0].set(ylabel='population 1', ylim=(0, subN))

   axs[1].plot(SME2.t / ms, SME2.i, '.', markersize=2, color='darkblue')
   axs[1].set(ylabel='population 2', ylim=(0, subN))

   axs[2].plot(R1.t / ms, R1.smooth_rate(window='flat', width=100 * ms) / Hz, color='darkred')
   axs[2].plot(R2.t / ms, R2.smooth_rate(window='flat', width=100 * ms) / Hz, color='darkblue')
   axs[2].set(ylabel='Firing rate (Hz)')

   axs[3].plot(E1.t / ms, E1.rates[0] / Hz, color='darkred')
   axs[3].plot(E2.t / ms, E2.rates[0] / Hz, color='darkblue')
   axs[3].set(ylabel='Input (Hz)', xlabel='Time (ms)')

   fig.align_ylabels(axs)
   fig.suptitle(f'coherence of random dots : {coh}', fontsize=16)
   plt.savefig(f"Q3_{cohid}")
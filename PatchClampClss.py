import os 
import efel
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
from ipfx.feature_extractor import (SpikeFeatureExtractor, SpikeTrainFeatureExtractor)
import pyabf

class PatchClamp():

    def __init__(self, DPath , T0, TF, S_Path):
        self.Path = DPath 
        self.SPath = S_Path
        self.DataABF = pyabf.ABF(DPath)   
        self.T0 = T0
        self.TF = TF 

#----------------------------------------------------------------------------------
    def SFX(self, **kwargs):
        DataF = []
        for x in self.DataABF.sweepList:
            self.DataABF.setSweep(x)
            T = self.DataABF.sweepX
            V = self.DataABF.sweepY
            C = self.DataABF.sweepC
            T0, TF = self.T0, self.TF

            SFX = SpikeFeatureExtractor(T0, TF,
                                        filter=None,
                                        dv_cutoff=20.0,
                                        thresh_frac=0.05,
                                        min_peak=0)
            SFX_Results = SFX.process(T, V, C)
            DataF.append(SFX_Results)
        Tab = pd.concat(DataF)

        plt.figure()
        plt.xlabel('Tiempo (s)')
        plt.ylabel('Voltaje (mV)')
        for x in self.DataABF.sweepList:
            self.DataABF.setSweep(x)
            plt.plot(self.DataABF.sweepX, self.DataABF.sweepY, alpha = 0.5, label = 'Sweep %d' %(x))
        plt.xlim(self.T0, self.TF)

        
        if kwargs['SwpNumber'] != None :
            SN = kwargs['SwpNumber']
            self.DataABF.setSweep(SN)
            plt.plot(self.DataABF.sweepX, self.DataABF.sweepY, linewidth =1, color = 'black')
        else:
            pass

        if kwargs['Save'] == 'Y':
            NameFig = input('Escriba el nombre con el que se guardará la gráfica: ')
            NameTab = NameFig

            NameFig = NameFig + '.png'
            NameTab = NameTab + '.csv'

            SFpath = os.path.join(self.SPath, NameFig)
            STPath = os.path.join(self.SPath, NameTab)
            Tab.to_csv(STPath)
            plt.savefig(SFpath)
        else:
            pass

        plt.show()
        print(Tab)

#-------------------------------------------------------------------------------------------------
    def STFX(self, **kwargs):
        if kwargs['NTraces'] == 'One' :
            self.DataABF.setSweep(sweepNumber=10, channel= 0)
            T = self.DataABF.sweepX
            V = self.DataABF.sweepY
            C = self.DataABF.sweepC
            currents = []

            T1 = int(400*self.DataABF.dataPointsPerMs)
            T2 = int(500*self.DataABF.dataPointsPerMs)
            Cu_Mean = np.average(self.DataABF.sweepC[T1:T2])

            T0, TF = self.T0, self.TF
            SFX = SpikeFeatureExtractor(start=T0, end=TF,
                                    filter=None,
                                    dv_cutoff=20.0,
                                    min_peak=0)
            SFX_Res = SFX.process(T, V, C)
            STFX = SpikeTrainFeatureExtractor(T0, TF)
            STFX_Res = STFX.process(T, V, C, SFX_Res)

            Tab = pd.DataFrame()
            Len = len(Tab)

            Tab.loc[Len, 'C_Step'] = Cu_Mean
            Tab.loc[Len, 'Adaptation'] = STFX_Res['adapt']
            Tab.loc[Len, 'Latency'] = STFX_Res['latency']
            Tab.loc[Len, 'ISI_CV'] = STFX_Res['isi_cv']
            Tab.loc[Len, 'ISI_Mean'] = STFX_Res['mean_isi']
            Tab.loc[Len, 'ISI_First'] = STFX_Res['first_isi']
            Tab.loc[Len, 'Firing_rate_Hz'] = STFX_Res['avg_rate']

            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax1.plot(self.DataABF.sweepX, self.DataABF.sweepY)
            ax1.plot(SFX_Res['peak_t'], SFX_Res['peak_v'], 'r.')
            ax1.plot(SFX_Res['threshold_t'], SFX_Res['threshold_v'], 'k.')
            ax1.set_ylabel('Voltage (mV)')
            ax2 = fig.add_subplot(212, sharex = ax1)
            ax2.plot(self.DataABF.sweepX, self.DataABF.sweepY, color = 'r')
            ax2.set_ylabel('Current (pA)')
            ax2.set_xlabel('Time (s)')
            ax2.axes.set_xlim(0,1)

            if kwargs['SwpNumber'] != None :
                SN = kwargs['SwpNumber']
                self.DataABF.setSweep(SN)
                plt.plot(self.DataABF.sweepX, self.DataABF.sweepY, linewidth =1, color = 'black')
            else:
                pass    

            if kwargs['Save'] == 'Y':
                NameFig = input('Escriba el nombre con el que se guardará la gráfica: ')
                NameTab = NameFig

                NameFig = NameFig + '.png'
                NameTab = NameTab + '.csv'

                SFpath = os.path.join(self.SPath, NameFig)
                STPath = os.path.join(self.SPath, NameTab)
                Tab.to_csv(STPath)
                fig.savefig(SFpath)
            else:
                pass

            plt.show()
            print(Tab)

        elif kwargs['NTraces'] == 'Sev' :

            Tab = pd.DataFrame()
            for n in self.DataABF.sweepList:
                self.DataABF.setSweep(n)
                T = self.DataABF.sweepX
                V = self.DataABF.sweepY
                C = self.DataABF.sweepC

                Currents =[]
                T1 = int(400*self.DataABF.dataPointsPerMs)
                T2 = int(500*self.DataABF.dataPointsPerMs)
                Cu_Mean = np.average(self.DataABF.sweepC[T1:T2])

                T0, TF = self.T0, self .TF

                SFX = SpikeFeatureExtractor(start=T0, end=TF, filter=None)
                SFX_Res = SFX.process(T, V, C)
                STFX = SpikeTrainFeatureExtractor(start=T0, end=TF)
                STFX_Res = STFX.process(T, V, C, SFX_Res)

                Len = len(Tab)
                Tab.loc[Len, 'Paso_Act'] = Cu_Mean
                Tab.loc[Len, 'Avg_Rate'] = STFX_Res['avg_rate']
                if STFX_Res['avg_rate'] > 0:
                    Tab.loc[Len, 'Adaptación'] = STFX_Res['adapt']
                    Tab.loc[Len, 'Latencia_s'] = STFX_Res['latency']
                    Tab.loc[Len, 'ISI_CV'] = STFX_Res['isi_cv']
                    Tab.loc[Len, 'ISI_Prom_ms'] = STFX_Res['mean_isi']*1000
                    Tab.loc[Len, 'ISI_First_ms'] = STFX_Res['first_isi']*1000
                
            Fig = plt.figure(figsize= (10,5))

            ax1 = Fig.add_subplot(1,2,1)
            Currents = []
            for i in self.DataABF.sweepList:
                self.DataABF.setSweep(i)
                Currents.append(np.average(self.DataABF.sweepC[T1:T2]))
            ax1.scatter(Currents, Tab.loc[:,'Avg_Rate'])
            ax1.set_xlabel('Current step (pA)')
            ax1.set_ylabel('Action Potentials')

            ax2 = Fig.add_subplot(2, 2, 2)
            for i in self.DataABF.sweepList:
                self.DataABF.setSweep(i)
                ax2.plot(self.DataABF.sweepX, self.DataABF.sweepY, alpha = 0.5, label = 'Sweep %d' % (i))
            ax2.set_ylabel('Membrane Voltage (mV)')
            

            ax3= Fig.add_subplot(2,2,4,sharex =ax2)
            for j in self.DataABF.sweepList:
                self.DataABF.setSweep(j)
                Current = self.DataABF.sweepC
                ax3.plot(self.DataABF.sweepX, Current)
            ax3.set_ylabel('Current (pA)')
            ax3.set_xlabel('Time (s)')

            if kwargs['SwpNumber'] != None :
                SN = kwargs['SwpNumber']
                self.DataABF.setSweep(SN)
                plt.plot(self.DataABF.sweepX, self.DataABF.sweepY, linewidth =1, color = 'black')
            else:
                pass

            if kwargs['Save'] == 'Y':
                NameFig = input('Escriba el nombre con el que se guardará la gráfica: ')
                NameTab = NameFig

                NameFig = NameFig + '.png'
                NameTab = NameTab + '.csv'

                SFpath = os.path.join(self.SPath, NameFig)
                STPath = os.path.join(self.SPath, NameTab)
                Tab.to_csv(STPath)
                fig.savefig(SFpath, dpi = 300)

            plt.show()
            print(Tab)

#----------------------------------------------------------------------
    def STF_Spreads(self, CSVPath):

        Data = np.loadtxt(fname= CSVPath, delimiter=',')
        Tab = pd.DataFrame()

        T = Data[:,0]/1000
        V = Data[:, 1:14]
        C = Data[:, 14:28]

        for i, v in enumerate(v.T):
            Current = C.T[i]
            T0, TF = self.T0, self.TF
            STFX = SpikeTrainFeatureExtractor(start= T0, end = TF)
            SFX = SpikeFeatureExtractor(start= T0, end= TF, filter= None)
            SpikesDF = SFX.process(t = T, v = V, i = C)
            STFX = STFX.process(t = T, v = V, i = Current, spikes_df=SpikesDF)

            #Tab
            Len = len(Tab)
            Tab.loc[Len, 'avg_rate'] = STFX['avg_rate']
            if STFX['avg_rate'] > 0:
                Tab.loc[Len, 'Adaptación'] = STFX['adapt']
                Tab.loc[Len, 'Latency'] = STFX['latency']
                Tab.loc[Len, 'ISI_CV'] = STFX['isi_cv']
                Tab.loc[Len, 'ISI_mean'] = STFX['mean_isi']
                Tab.loc[Len, 'ISI_first'] = STFX['first_isi']

            plt.plot(T, V)
            plt.xlim(0, 1)
            plt.ylabel('Tiempo (s)')
            plt.xlabel('Voltaje (mV)')
        
        t1 = C[(T > 0.3) & (T < 0.4), :]
        t2 = C[(T > 0) & (T < 0.2), :]
        C_Mean = np.median(t1, axis = 0) -np.median(t2, axis = 0)
        Tab['Current_steps'] = np.round(C_Mean, 0)
        plt.show()
        print(Tab)
#------------------------------------------------------------------------------------------------

    def SFX_Efel(self, SN):
        Sweep = self.DataABF.setSweep(sweepNumber=SN, channel=0)

        T = self.DataABF.sweepX*1000
        V = self.DataABF.sweepY
        C = self.DataABF.sweepC

        Trace = {'T': self.DataABF.sweepX*1000,
                 'V': self.DataABF.sweepY,
                 'stim_start': [100],
                 'stim_end': [800]} 
        Traces = [Trace]

        Currents = []
        T1 = int(400*self.DataABF.dataPointsPerMs)
        T2 = int(500*self.DataABF.dataPointsPerMs)
        C_Mean = np.average(self.DataABF.sweepC[T1, T2])

        efel.api.setThreshold(0)
        efel.api.setDerivativeThreshold(20)

        Feat_Values = efel.getFeatureValues(Traces, 
                                            ['AP_amplitude', 'AP_witdh',
                                             'AP_begin_voltage','AP_rise_rate',
                                             'AP_fall_rate'], raise_warnings=None)[0]
        
        Tab = pd.DataFrame(columns=['Sweep', 'Current_step', 
                            'AP_amplitude', 'AP_width', 
                            'AP_begin_voltage', 'AP_rise_rate', 'AP_fall_rate'])
        
        Len = len(Tab)
        Tab.loc[Len, 'Sweep'] = SN
        Tab.loc[Len, 'Current_Step'] = C_Mean
        Tab.loc[Len, 'AP_Amplitude'] = Feat_Values['AP_amplitude'][0]
        Tab.loc[Len, 'AP_width'] = Feat_Values['AP_width'][0]
        Tab.loc[Len, 'AP_begin_voltage'] = Feat_Values['AP_begin_voltage'][0]
        Tab.loc[Len, 'AP_rise_rate'] = Feat_Values['AP_rise'][0]
        Tab.loc[Len, 'AP_fall_rate'] = Feat_Values['AP_fall_rate'][0]

        Fig = plt.figure()
        ax1 = Fig.add_subplot(121)
        ax1.plot(T, V)
        ax1.set_xlabel('Tiempo (ms)')
        ax1.set_ylabel('Voltaje de membrana (mV)')
        ax1.axes.set_xlim(0,1000)

        ax2 = Fig.add_subplot(122, sharey =ax1)
        ax2.plot(T, V)
        ax2.set_xlabel('Tiempo (ms)')
        ax2.set_ylabel('Voltaje de membrana (mV)')
        ax2.axes.set_xlim(240, 254)
        plt.tight_layout()
        plt.show()

        plt.show()
        print(Tab)
        
#-------------------------------------------------------------------------------------------------    
    def STFX_Efel(self, SN):
        
        Tab = pd.DataFrame(columns=['Current_Step', 'Spikecount', 'adaptation',
                                    'Latency_ms', 'ISI_CV', 'ISI_mena_ms'])
        
        for i in self.DataABF.sweepList:
            self.DataABF.setSweep(i)

            T = self.DataABF.sweepX*1000
            V = self.DataABF.sweepY

            Trace = {'T': self.DataABF.sweepX*1000,
                     'V': self.DataABF.sweepY,
                     'stim_start': [100], 'stim_end': [800]}
            Traces = [Trace]

            efel.api.setThreshold(0)
            efel.api.setDerivativeThreshold(20)

            Feat_Val = efel.getFeatureValues(Traces, 
                                    ['Current_Step', 'Spikecount', 'adaptation',
                                    'Latency_ms', 'ISI_CV', 'ISI_mena_ms'],
                                     raise_warnings=None)[0]
            
            C =self.DataABF.sweepC
            Cs = []
            T0 = int(400*self.DataABF.dataPointsPerMs)
            TF = int(500*self.DataABF.dataPointsPerMs)
            C_Mean = np.average(self.DataABF.sweepC[T0:TF])
            Cs.append(C_Mean)

            Len = len(Tab)
            Tab.loc[Len, 'Current_step'] = C_Mean
            Tab.loc[Len, 'Spikecount'] = Feat_Val['Spikecount'][0]
            if Feat_Val['Spikecount'] is not None:
                Tab.loc[Len, 'Latency_ms'] = Feat_Val['time_to_first_spike']
                if Feat_Val['Spikecount'] > 4:
                    Tab.loc[Len, 'adaptation'] = Feat_Val['adaptation_index'][0]
                    Tab.loc[Len, 'ISI_CV'] = Feat_Val['ISI_CV'][0]
                    Tab.loc[Len, 'ISI_mean'] = Feat_Val['Isi_Values'][0]/1000

        Fig = plt.figure()
        ax1= Fig.add_subplot(121)
        for i in self.DataABF.sweepList:
            self.DataABF.setSweep(i)
            ax1.plot(self.DataABF.sweepX*1000, self.DataABF.sweepY, alpha=0.3)
        ax1.set_xlabel('Tiempo (ms)')
        ax1.set_ylabel('Membrane voltaje (mV)')
        ax1.axes.set_xlim(0, 1000)

        if SN != None:
            ax2 = Fig.add_subplot(122, sharey = ax1)
            ax2.plot(T, V, label = 'Sweep %d' %(SN))
            ax2.set_xlabel('Tiempo (ms)')
            ax2.set_ylabel('Membrane voltaje (mV)')
            ax2.set_xlim(0,1000)
            plt.tight_layout()
            plt.show()
            print(Tab)
        else: 
            plt.show() 
            print(Tab)     
#--------------------------------------------------------------------------------------------------

    def STFX_NA(self, TH, Save_T, F_Name):
        Data = np.loadtxt(fname=self.Path, delimiter=',')

        T = Data[:0]
        V = Data[:,1:14]
        C = Data[:, 14:28]

        Tab = pd.DataFrame(columns= ['Spikecount', 'adaptation','Latency_ms',
                                     'ISI_CV', 'ISI_mean_ms'])
        
        for v in V.T:
            Trace ={'T':T, 'V': V, 'stim_start' :[100],
                    'steam_end': [700]}
            Traces = [Trace]
            efel.api.setThreshold(0)
            efel.api.setDerivativeThreshold(20)

            Feat_V = efel.getFeatureValues(Traces,['Spikecount', 'adaptation_index', 'time_to_first_spike', 
                                                   'ISI_CV', 'ISI_Values'], raise_warnings=None)[0]
            Len = len(Tab)
            Tab.loc[Len, 'Spikecount'] = Feat_V['Spikecount'][0]

            if Feat_V['Spikecount'] is not None:
                Tab.loc[Len, 'Latency_ms'] = Feat_V['time_to_first_spike']
                if Feat_V['Spikecount'] >4:
                    Tab.loc[Len, 'adaptation'] = Feat_V['adaptation_index'][0]
                    Tab.loc[Len, 'ISI_cv'] = Feat_V['ISI_CV'][0]
                    Tab.loc[Len, 'ISI_mean_ms'] = Feat_V['ISI_Values'][0]/1000

        T0 = C[(T > 300) & (T<400 ), :]
        TF = C[(T > 0) & (T < 200 ), :]
        C_Mean = np.median(T0, axis = 0) - np.median(TF, axis = 0)
        Tab['Current_steps'] = np.round(C_Mean, 0)

        Fig = plt.figure()
        plt.plot(T, V, linewidth = 1, alpha = 0.6)

        if TH != None :
            Trace_H = Data[:,4]
            plt.plot(T, Trace_H, linewiidth = 1, color= 'b' )
        else:
            pass
        
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")
        plt.xlim(0, 1000)
        plt.show()
        print(Tab)

        if Save_T == 'Y' & F_Name != None :
            C_WS = os.path.join(self.SPath, F_Name + '.csv')
            Tab.to_csv(C_WS)

#----------------------------------------------------------
    def PhaseS(self, Data, Swp, SR, HW_ms, Save):
        ABFData = pyabf.ABF(Data)

        def dvdt(abf, Sweep, sampling_rate, half_window_ms):

            sweep = abf.setSweep(sweepNumber =Sweep, channel = 0)
            dv = np.diff(abf.sweepY)
            dv_Max = np.argmax(dv)
            SR = sampling_rate
            HWms = half_window_ms
            HW = (HWms*SR)//1000
            global t_window
            t_window = abf.sweepX[dv_Max-HW : dv_Max+ HW]
            v_window =  abf.sweepY[dv_Max-HW : dv_Max+ HW]
            dv_ap = np.diff(v_window)
            dt_ap = np.diff(t_window*1000)
            dvdt = dv_ap/dt_ap
            dvdtmax = np.amax(dvdt)
            v_ap_array = np.delete(v_window, 1)
            return {'voltage': v_window, 'dvdt': dvdt,
                    'dvdtMax':dvdtmax, 'time':t_window,
                    'voltage_array': v_ap_array}
        
        dvdt = dvdt(ABFData, Swp, SR, HW_ms)
        Tab = pd.DataFrame(columns=['voltage','dvdt'])
        Tab.voltage = pd.DataFrame(dvdt['voltage'])
        Tab.dvdt = pd.DataFrame(dvdt['dvdt'])

        Fig = plt.figure()
        ax1 = Fig.add_subplot()
        ax1.set_title('Phase-plane plot')
        ax1.plot(dvdt['time'], dvdt['voltage'])
        ax1.set_ylabel('V (mV)')
        ax1.set_xlabel('t (s)')

        ax2 = Fig.add_subplot(122)
        ax2.set_title('Phase-plane plot')
        ax2.set_ylabel('dv/dt (mV/ms)')
        ax2.set_xlabel('V (mV)')

        print('dv/dt Max =', dvdt['dvdtMax'])
        Fig.tight_layout()
        plt.show()
        print(Tab)
        
        if Save == 'Y':

            FileN = input('Introduzca el nombre con el')

            Fig.savefig(self.SPath + FileN + '.png', dpi = 300)
            Tab.to_csv(self.SPath + FileN + '.csv')
        
        else: 
            pass 


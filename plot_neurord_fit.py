import sys
ajupath='C:\\Users\\kblackw1\\Documents\\Python_Scripts\\ajustador'
if ajupath not in sys.path:
    sys.path.insert(0,ajupath)
import ajustador as aju
from ajustador import loadconc,nrd_fitness,nrd_output
import numpy as np
import os
import glob
from matplotlib import pyplot as plt
plt.ion()

#FIXME: put loadnpz and plot_history into utilities, with additional_info passed in
#then, import those functions into both plot_opt_waves.py and plot_neurord_fit.py
def load_npz(fname,additional_info):
    dat = np.load(fname)
    fitvals=dat['fitvals']
    tmpdirs=dat['tmpdirs']
    paramvals=dat['params']
    paramnames=dat['paramnames']
    features=dat['features']
    info_dict={}
    for feat in additional_info:
        print(dat['features'])
        text=[item for item in dat['features'] if item.startswith(feat+'=')][0]
        if len(text):
            info_dict[feat]=text.split(feat+'=')[-1]
    return fitvals,tmpdirs,paramvals,paramnames,info_dict

def plot_history(fitvals,dataname):
    plt.figure()
    plt.plot(fitvals[:,-1],'r.')
    plt.title(dataname)
    plt.xlabel('evaluation')
    plt.ylabel('fitness')
    return

def plot_traces(bestdir,dataname,traces,exp_to_fit):
    colors=[plt.get_cmap('viridis'),plt.get_cmap('gist_heat')]#,plt.get_cmap('summer')]
    color_increment=int(round(256.0/(len(exp_to_fit.data)+1)))
    base=int(color_increment/2)
    fig,axes=plt.subplots(nrows=len(traces[0].keys()),ncols=1)
    for dd in range(len(exp_to_fit.data)):
        dcolor = colors[0].__call__(base+dd*color_increment % colors[0].N)
        scolor = colors[1].__call__(base+dd*color_increment % colors[0].N)
        print(dcolor,scolor)
        for ax,(mol,csvconc) in enumerate(exp_to_fit.data[dd].waves.items()):
            axes[ax].plot(csvconc.wave.x,csvconc.wave.y,label='data '+str(dd),color=dcolor)
            axes[ax].plot(traces[dd][mol]['x'],traces[dd][mol]['y'],label='sim '+str(dd),color=scolor)
            axes[ax].set_ylabel(mol+' Conc (nM)')
        axes[ax].legend()
    fig.suptitle(dataname)
    axes[ax].set_xlabel('Time (msec)')
    return

if __name__ == '__main__':
    #directories with npz file and tmpdirs
    npzfit_dir='/home/kblackw1/sigpath/cofilin/cof_fit/cof_opt_sum/'
    #'C:\\Users\\kblackw1\\Documents\\HippoAkap\\cofilin\\cofilin\\cof_opt\\'
    #newdirname='C:\\Users\\kblackw1\\Documents\\HippoAkap\\cofilin\\cofilin\\ab_cof_opt\\'
    #parameters from fit.py file
    exp_name='Bosch_Hedrick_cof-basal300'
    mol={'RacGTP':['RacGTP'],'Cofactin':['Cofactin','Cof']}
    start_stim=90
    #norm='percent'
    # can have 'neuron' if fitting to moose_nerp
    additional_info=['model','dataname']

    ### end parameter specification
    os.chdir(npzfit_dir)
    pattern='*.npz'
    fnames=glob.glob(pattern)

    if len(fnames)==1:
        fitvals,tmpdirs,params,paramnames,mod_neur_dict=load_npz(fnames[0],additional_info)
        exp_data=loadconc.CSV_conc_set(exp_name, stim_time=start_stim)
        plot_history(fitvals,exp_name)
        bestfit=np.argmin(fitvals[:,-1])
        while bestfit>-1:
            bestdir=tmpdirs[bestfit] #Edit if using NSG - see plot_opt_results.py
            if not os.path.exists(bestdir): #update tmpdirs with new directory if fit dirs have been moved/copied
                dirname=os.path.basename(bestdir)
                bestdir=newdirname+dirname
                #TEMPORARY
                #bestdir=newdirname+'tmp8o0z7q6l\\'
            bestfitval=fitvals[bestfit,-1]
            simfiles=glob.glob(bestdir+'/*.h5')
            if len(simfiles)>0:
                print('h5 file=',os.path.basename(simfiles[0]), ', molecules=', exp_data.data[0].waves.keys(),',bestfit number=',bestfit)
                traces=[{sp:{} for sp in mol.keys()} for f in simfiles]
                for dd,f in enumerate(simfiles):
                    result=nrd_output.Output(f,start_stim)
                    for sp in mol.keys():
                        if exp_data.data[0].waves[sp].norm:
                            traces[dd][sp]['y'],traces[dd][sp]['x']=nrd_fitness.nrd_output_percent(result,mol[sp],start_stim,exp_data.data[0].waves[sp])
                        else:
                            traces[dd][sp]['y'],traces[dd][sp]['x']=nrd_fitness.summed_species(result,mol[sp])
                plot_traces(bestdir,exp_name,traces,exp_data)
                previous_bestfit=bestfit
                text=input('do you want to try another fit (fitnumber for yes /-1 for no)')
                bestfit=int(text)
            else:
                print('No h5 file found in ', bestdir)
        bestfit=previous_bestfit
        print ('model parameters for', previous_bestfit, 'in ',bestdir+'model-0.xml')
            

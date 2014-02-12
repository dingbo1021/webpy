import web
import math
import os
import re
import random
import Image
import matplotlib.pyplot as plt1
import matplotlib.pyplot as plt2
import matplotlib.pyplot as plt3
from web import form
render = web.template.render('templates')

rootpath = 'C:/webpy/'

urls = (
    '/startpage', 'startpage',
    '/customize', 'customize',
    '/para', 'para',
    '/result', 'result',
    '/result/(.*)','checkback'
 #   '/result/\d{11,13}','checkback'
)
app = web.application(urls, globals())

myform = form.Form( 
    form.Dropdown('Pfile selection:', 
    [('HyMap_Test_Web_Tree.fcm','Predefined scenario BKG TREE'), 
    ('Hyp_Web_Test.fcm','Predefined scenario BKG GRASS'), 
    ('HyMap_Test_Web.fcm','Predefined scenario 3 BKGs'),
    ('Customized.fcm','User Customize')]
    )) 

class startpage:
    def GET(self):
        form = myform()
        return render.startpage(form)
    def POST(self): 
        form = myform() 
        if not form.validates(): 
            return render.startpage(form)
        else:
            if form['Pfile selection:'].value == 'HyMap_Test_Web_Tree.fcm':
                return render.case1()
            if form['Pfile selection:'].value == 'Hyp_Web_Test.fcm':
                return render.case2()
            if form['Pfile selection:'].value == 'HyMap_Test_Web.fcm':
                return render.case3()
            if form['Pfile selection:'].value == 'Customized.fcm':
 #               raise web.seeother('/customize')
                return render.customize()
                
class customize:
    def GET(self):
        return render.jump()
        
    def POST(self): 
        bkg = web.input().get('hiddentxt') 
        sensor = web.input().get('sensor_selection') 
        target = web.input().get('target_selection') 

        bkg = bkg[:-1].strip().split(',')

        bkgall=''
        WVall = ''
        WV = readWVfile(sensor)
        percentall='100'+(len(bkg)-1)*',0'

#generate prefill background names and wavelength list
        for i in range(0, len(bkg)):
            bkgall = bkgall + '<option value="'+str(i)+'">'+bkg[i]+'</option>'
        for i in range(0, len(WV)):
            WVall = WVall + '<option value="'+str(i)+'">'+str(WV[i])+'</option>'

        return render.para(bkgall,','.join(bkg),sensor,target,WVall,percentall)

class para:
    def GET(self):
        return render.jump()
        
    def POST(self): 
        hiddentxt = web.input().get('hiddentxt') #the chosen wavelength
#        raise Exception(web.input())     
        Atmospheric_haze = web.input().get('Atmospheric_haze')
        Solangle = web.input().get('Solangle') 
        Atmospheric_model = web.input().get('Atmospheric_model')
        CloadIndex = web.input().get('icld')
        Ground_altitude = web.input().get('Ground_altitude')
        Meteorological_Range = web.input().get('Mrange')
        
        SensorName = web.input().get('Sensorname')
        Noisefac = web.input().get('Noisefac') 
        Gainfac = web.input().get('Gainfac') 
        Relcal = web.input().get('Relcal') 
        Platalt = web.input().get('Platalt') 
        Tint = web.input().get('Tint') 
        
        Targscale = web.input().get('Targscale') 
        Targperc = web.input().get('Targperc') 
        Targinback = web.input().get('Targinback')
        Targinback = str(int(Targinback)+1)
        TargName = web.input().get('Targname')
   
        Bkgscale = web.input().get('Bkgscale') 
        Backperc = web.input().get('Backperc') 
        Bkgname = web.input().get('Bkgname2') 
      
        while not Backperc[-1].isdigit():
            Backperc = Backperc[:-1]

#create a sessionkey for diff session
        session_key = str(random.random())[2:]
        sessionpath = rootpath+'static/'+session_key+'/'
        if os.path.exists(sessionpath)==False:
            os.mkdir(sessionpath);  

#modify custimized fcm file
        FCMdata = open(rootpath+'resource/pfiles/customized_base.fcm').read()

        FCMdata = modifyCustomFCM('ihaze', str(Atmospheric_haze), FCMdata)
        FCMdata = modifyCustomFCM('solangle', str(Solangle), FCMdata)
        FCMdata = modifyCustomFCM('gndalt', str(Ground_altitude), FCMdata) 
        FCMdata = modifyCustomFCM('model', str(Atmospheric_model), FCMdata) 
        FCMdata = modifyCustomFCM('icld', str(CloadIndex), FCMdata)  
        FCMdata = modifyCustomFCM('metrange', str(Meteorological_Range), FCMdata)              

        FCMdata = modifyCustomFCM('sensorfile', str(SensorName), FCMdata)
        FCMdata = modifyCustomFCM('noisefac', str(Noisefac), FCMdata)
        FCMdata = modifyCustomFCM('gainfac', str(Gainfac), FCMdata)
        FCMdata = modifyCustomFCM('relcal', str(Relcal), FCMdata)        
        FCMdata = modifyCustomFCM('platalt', str(Platalt), FCMdata)
        FCMdata = modifyCustomFCM('tint', str(Tint), FCMdata)

        FCMdata = modifyCustomFCM('targscale', str(Targscale), FCMdata)            
        FCMdata = modifyCustomFCM('targperc', str(Targperc), FCMdata,1)
        FCMdata = modifyCustomFCM('targinback', str(Targinback), FCMdata)
        FCMdata = modifyCustomFCM('targetname', str(TargName), FCMdata)

        FCMdata = modifyCustomFCM('backscale', str(Bkgscale), FCMdata)
        FCMdata = modifyCustomFCM('backname', str(Bkgname), FCMdata,1)
        FCMdata = modifyCustomFCM('backperc', Backperc, FCMdata,1)
        FCMdata = modifyCustomFCM('numback', str(len(Backperc.split(','))), FCMdata)

        ##########writing chosen wvlg to FCM file for C++ code
        FCMdata = modifyCustomFCM('wavelengthchosen', ','.join(str(hiddentxt).split('\r\n'))[:-1], FCMdata,1)

        open(rootpath+'resource/pfiles/customized.fcm', 'wb').write(FCMdata)
        open(rootpath+'static/'+session_key+'/customized.fcm', 'wb').write(FCMdata)

        ############################################################################################################

        feedback = os.system(rootpath+'working.exe')
        #to be finished: give fcm file to queue, then check sessionid folder. 
        ############################################################################################################   
#        feedback =0;   
        FCMdata=''
        if feedback == -1:
            return render.error()
        Wvlength = readfilep(str(SensorName),rootpath+'resource/sensorWV/')
        WvlengthCount = len(Wvlength)
        
#read txt file generated by working.exe
        Ltnumbers = readfilep('LBT')
        SNRn = readfilep('SNR')
        ROC = readfilep('ROC')

#generate checkbox for current background and target 
        pltselectR = ''
        pltselectS = ''
        for i in range(0, len(Bkgname.split(','))):
            pltselectR = pltselectR + '<li><label for="id_BkgnameR_'+str(i)+'"><input type="checkbox" name="BkgnameR" value="'+str(i)+'" id="id_BkgnameR_'+str(i)+'" /> '+Bkgname.split(',')[i]+'</label></li>'
            pltselectS = pltselectS + '<li><label for="id_BkgnameS_'+str(i)+'"><input type="checkbox" name="BkgnameS" value="'+str(i)+'" id="id_BkgnameS_'+str(i)+'" /> '+Bkgname.split(',')[i]+'</label></li>'
            
        i=i+1;
        pltselectR = pltselectR + '<li><label for="id_BkgnameR_'+str(i)+'"><input type="checkbox" name="BkgnameR" value="'+str(i)+'" id="id_BkgnameR_'+str(i)+'" /> '+'Avg Background'+'</label></li>'
        pltselectS = pltselectS + '<li><label for="id_BkgnameS_'+str(i)+'"><input type="checkbox" name="BkgnameS" value="'+str(i)+'" id="id_BkgnameS_'+str(i)+'" /> '+'Avg Background'+'</label></li>'
            
        i=i+1;
        pltselectR = pltselectR + '<li><label for="id_BkgnameR_'+str(i)+'"><input type="checkbox" name="BkgnameR" value="'+str(i)+'" id="id_BkgnameR_'+str(i)+'" /> Target: '+str(TargName)+'</label></li>'
        pltselectS = pltselectS + '<li><label for="id_BkgnameS_'+str(i)+'"><input type="checkbox" name="BkgnameS" value="'+str(i)+'" id="id_BkgnameS_'+str(i)+'" /> Target: '+str(TargName)+'</label></li>'

#plot radiance, SNR, ROC, default is to plot target radiance and SNR
        TgtName = []
        TgtName.append('Target: '+str(TargName))
        showback = len(Ltnumbers)/WvlengthCount-1 # default, target radiance and SNR are the last in data file
        
        LTdisplay = Ltnumbers[(showback*WvlengthCount):(showback*WvlengthCount+WvlengthCount)]    
        plt1.plot(Wvlength,LTdisplay)
        plt1.xlim([400, 2500]) 
        plt1.title('Scene Mean Spectral Radiance')
        plt1.legend(TgtName,prop={'size':8})
        plt1.xlabel('Wavelength(microns)')
        plt1.ylabel('Spectral Radiance(mW/cm^2-sr-um)')
        plt1.savefig(rootpath+'static/'+session_key+'/rad.png',dpi=100)
        plt1.clf()

        SNRdisplay = SNRn[showback*WvlengthCount:(showback*WvlengthCount+WvlengthCount)]
        plt2.plot(Wvlength,SNRdisplay)
        plt2.xlim([400, 2500]) 
        plt2.ylim([0, 110]) 
        plt2.title('Sensor Signal-to-Noise Ratio')
        plt2.legend(TgtName,prop={'size':8})
        plt2.xlabel('Wavelength(microns)')
        plt2.ylabel('Signal-to-Noise Ratio')
        plt2.savefig(rootpath+'static/'+session_key+"/snr.png",dpi=100)
        plt2.clf()

        Pdmin = ROC[0:73]
        Pfa = ROC[73:146]
        
        plt3.semilogx(Pfa,Pdmin)            
        plt3.plot(Pfa,Pdmin,"b")
        plt3.axis([1e-4, 1, 0, 1.1])
        plt3.title('ROC Curve')
        plt3.xlabel('Probability of False Alarm')
        plt3.ylabel('Probability of Detection')
        plt3.savefig(rootpath+'static/'+session_key+'/roc.png',dpi=100)
        plt3.clf()   
        
        
        f = open(rootpath+'LBT.txt', 'r')
        dataL = f.readlines()
        f.close()
        
        f = open(rootpath+'SNR.txt', 'r')
        dataSNR = f.readlines()
        f.close()
        
        f = open(rootpath+'ROC.txt', 'r')
        dataROC = f.readlines()
        f.close()
        
        DataTitleL = [];
        DataTitleSNR = [];
        
#generate radiance, SNR, ROC with title line
        datasize = len(dataL);
        for x in range(0,datasize/WvlengthCount-2):
            DataTitleL.append('L of Background ' + str(x+1) + '\n');
            DataTitleSNR.append('SNR of Background ' + str(x+1) + '\n');
        DataTitleL.append('L of Background Average\n');
        DataTitleL.append('L of Target\n');
        DataTitleSNR.append('SNR of Background Average\n');
        DataTitleSNR.append('SNR of Target\n');

        for x in range(0,datasize/WvlengthCount):
            dataL.insert(datasize-WvlengthCount*(x+1),DataTitleL[datasize/WvlengthCount-1-x])
            dataSNR.insert(datasize-WvlengthCount*(x+1),DataTitleSNR[datasize/WvlengthCount-1-x])

        dataROC.insert(0,'P detection\n')
        dataROC.insert(74,'P false alarm\n')
        
#write LBT, SNR, ROC file to corresponding folder
        f = open(sessionpath+'ROC.txt', 'wb')
        f.write(''.join(dataROC))
        f.close()
        f = open(sessionpath+'LBT.txt', 'wb')
        f.write(''.join(dataL))
        f.close()
        f = open(sessionpath+'SNR.txt', 'wb')
        f.write(''.join(dataSNR))
        f.close()

        return render.result(SensorName,session_key,pltselectR,pltselectS,str(TargName),str(Bkgname))
class result:
    def GET(self):
        return render.jump()
        
    def POST(self): 
        

        Radiance=[]
        SNR=[]
        choices=[]
        
        session_key = web.input().get('session_key')
        SensorName = web.input().get('SensorName')
        Bkgname = web.input().get('bkgname').split(',')
        Tgtname = web.input().get('tgtname')

#get corresponding data file according to sensor name and session
        Wvlength = readfilep(SensorName, rootpath+'resource/sensorWV/')
        WvlengthCount = len(Wvlength)
        
        sessionpath = rootpath+'static/'+session_key+'/'
        Ltnumbers = readfilepstr('LBT',sessionpath)
        SNRn = readfilepstr('SNR',sessionpath)
        ROC = readfilep('ROC')
        
#generate checkbox for current background and target 
        pltselectR = ''
        pltselectS = ''

        for i in range(0, len(Bkgname)):
            pltselectR = pltselectR + '<li><label for="id_BkgnameR_'+str(i)+'"><input type="checkbox" name="BkgnameR" value="'+str(i)+'" id="id_BkgnameR_'+str(i)+'" /> '+Bkgname[i]+'</label></li>'
            pltselectS = pltselectS + '<li><label for="id_BkgnameS_'+str(i)+'"><input type="checkbox" name="BkgnameS" value="'+str(i)+'" id="id_BkgnameS_'+str(i)+'" /> '+Bkgname[i]+'</label></li>'
            
        i=i+1;
        pltselectR = pltselectR + '<li><label for="id_BkgnameR_'+str(i)+'"><input type="checkbox" name="BkgnameR" value="'+str(i)+'" id="id_BkgnameR_'+str(i)+'" /> '+'Avg Background'+'</label></li>'
        
        pltselectS = pltselectS + '<li><label for="id_BkgnameS_'+str(i)+'"><input type="checkbox" name="BkgnameS" value="'+str(i)+'" id="id_BkgnameS_'+str(i)+'" /> '+'Avg Background'+'</label></li>'
            
        i=i+1;
        pltselectR = pltselectR + '<li><label for="id_BkgnameR_'+str(i)+'"><input type="checkbox" name="BkgnameR" value="'+str(i)+'" id="id_BkgnameR_'+str(i)+'" /> Target: '+Tgtname+'</label></li>'
        
        pltselectS = pltselectS + '<li><label for="id_BkgnameS_'+str(i)+'"><input type="checkbox" name="BkgnameS" value="'+str(i)+'" id="id_BkgnameS_'+str(i)+'" /> Target: '+Tgtname+'</label></li>'
                 
#get checkbox values for radiance and SNR and parse to index number for redraw                 
        A=web.data().split('&')        
        A=A[2:-2]
        B='.'.join(A)
        B=B.replace('BkgnameR','R')
        B=B.replace('BkgnameS','S')

        while B.find('R')!=-1:
            Radiance.append(int(B[B.find('R')+2]))
            B=B[4:]
        while B.find('S')!=-1:
            SNR.append(int(B[B.find('S')+2]))
            B=B[4:]
            
        choices = Bkgname
        choices.append('Avg Background')
        choices.append('Target:'+Tgtname)
#redraw Radiance if necessary
        if len(Radiance) != 0:
            LegendRadiance=[]
            for i in range(0, len(Radiance)):
                LegendRadiance.append(choices[int(Radiance[i])])
                LTdisplay = Ltnumbers[(int(Radiance[i])*(WvlengthCount+1)+1):(int(Radiance[i])*(WvlengthCount+1)+1+WvlengthCount)]    
                plt1.plot(Wvlength,map(float,LTdisplay))
#            raise Exception(LegendRadiance)    
            plt1.xlim([400, 2500]) 
            plt1.title('Scene Mean Spectral Radiance')
            plt1.xlabel('Wavelength(microns)')
            plt1.ylabel('Spectral Radiance(mW/cm^2-sr-um)')
            plt1.legend(LegendRadiance,prop={'size':8})
            plt1.savefig(rootpath+'static/'+session_key+'/rad.png',dpi=100)
            plt1.clf()
#redraw SNR if necessary  
        if len(SNR) != 0:
            LegendSNR=[]  
            for i in range(0, len(SNR)):
                LegendSNR.append(choices[int(SNR[i])])
                SNRdisplay = SNRn[(int(SNR[i])*(WvlengthCount+1)+1):(int(SNR[i])*(WvlengthCount+1)+1+WvlengthCount)]
                plt2.plot(Wvlength,map(float,SNRdisplay))
            plt2.xlim([400, 2500]) 
            plt2.ylim([0, 105]) 
            plt2.title('Sensor Signal-to-Noise Ratio')
            plt2.xlabel('Wavelength(microns)')
            plt2.ylabel('Signal-to-Noise Ratio')
            plt2.legend(LegendSNR,prop={'size':8})
            plt2.savefig(rootpath+'static/'+session_key+'/snr.png',dpi=100)
            plt2.clf()           
                

        return render.result(SensorName,session_key,pltselectR,pltselectS, Tgtname, web.input().get('bkgname'))
class checkback:
    def GET(self,session_key):
        if os.path.exists(rootpath+'static/'+session_key+'/')==False:
            return web.notfound('4')
        else:
            return render.cbresult(session_key)
if __name__ == "__main__":
    app.run()

def readWVfile(filename):
    str1 = rootpath+'resource/sensorWV/'
    str2 = filename;
    str3 = '.txt'
    txt = open(str1+str2+str3,'r')
    txtall = txt.readlines()
    Output = map(float,txtall)    
    txt.close()
    return Output
def readfilep(filename,path=rootpath):
    str2 = filename;
    str3 = '.txt'
    txt = open(path+str2+str3,'r')
    txtall = txt.readlines()
    Output =map(float,txtall)
    txt.close()
    return Output
    
def modifyCustomFCM(keyword, value, data, multi = 0):
    if value.count(',')==0 and multi==0:
        data = re.sub(keyword + r'.*' +';', keyword + '\t' + value +'\t'+';',data)
    else:
        data = re.sub(keyword + r'.*' +';', keyword + '\t/' + value +'/\t'+';',data)
    return data
 
def readfilepstr(filename,path):
    str2 = filename;
    str3 = '.txt'
    txt = open(path+str2+str3,'r')
    txtall = txt.readlines()
    Output =map(str,txtall)    
    return Output    

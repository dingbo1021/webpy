#include "classes.h"

#define DEBUG

int main(int argc, char *argv[])
{ 
	try
	{
	////predefined pathconfig file path
    FileDir* m_FileDir = new FileDir("C:/webpy/pathconfig.ini");
//    char* dir = "HyMap_Test_Web.fcm";
    m_FileDir->showall();

	ofstream fLBT(m_FileDir->getLBT());  
	ofstream fSNR(m_FileDir->getSNR());  
	ofstream fROC(m_FileDir->getROC());
	ofstream fdebug(m_FileDir->getForDebug());

	if(!(fLBT.is_open()&&fSNR.is_open()&&fROC.is_open()&&fdebug.is_open()))
	{
		logfile("Output file extablish failure.","/ERROR");
		throw GeneralException();
	}

	//// \\or// Linux problem
    string pfilepath = m_FileDir->getpfile() + "/" + "customized.fcm";
    
    Background Bkg(pfilepath.c_str());
    Target Tgt(pfilepath.c_str());
    Scene Sce(pfilepath.c_str());
    Sensor Sen(pfilepath.c_str(),m_FileDir);
    Processing Pro(pfilepath.c_str());
    Modtran Modtran(m_FileDir);

	Sen.showall();


	vector<int> wchosen = Sen.getbackchosen();
	if(0==wchosen.size())
	{
		logfile("Wavelengths are not chosen","/ERROR");
		throw GeneralException();
	}
    sort(wchosen.begin(),wchosen.end());


	int bkgCount = Bkg.getnumback();
    
    int SensorDims = Sen.getWlgth().size();
    int FeatureDims = wchosen.size();
    //  int SourceDims = 210;   
    int extrarun = 4; 
    int runsize = bkgCount+extrarun;
    
    MatrixXd FeatureMtrx(FeatureDims, SensorDims);

    FeatureMtrx.setZero();
    for(int i=0;i<FeatureDims;i++)
        FeatureMtrx(i,wchosen[i]) = 1;
    
    ArrayXd* WVNo = new ArrayXd[runsize];

//	ArrayXd* Trans = new ArrayXd[runsize];
	ArrayXd* SurfEmis = new ArrayXd[runsize];
	ArrayXd* SolScat = new ArrayXd[runsize];
	ArrayXd* GndRflt = new ArrayXd[runsize];
//	ArrayXd* DrctRflt = new ArrayXd[runsize];
//	ArrayXd* TotRad = new ArrayXd[runsize];

	//ArrayXd* TransO = new ArrayXd[runsize];
	ArrayXd* SurfEmisO = new ArrayXd[runsize];
	ArrayXd* SolScatO = new ArrayXd[runsize];
	ArrayXd* GndRfltO = new ArrayXd[runsize];
	//ArrayXd* DrctRfltO = new ArrayXd[runsize];
	//ArrayXd* TotRadO = new ArrayXd[runsize];
    
	ArrayXd* MeansB = new ArrayXd[bkgCount+1];
    MatrixXd* CovsB = new MatrixXd[bkgCount+1];
    MatrixXd* CovsLB = new MatrixXd[bkgCount+1];
	ArrayXd* NoiseB = new ArrayXd[bkgCount+1];
	ArrayXd* LB = new ArrayXd[bkgCount+1];
	ArrayXd* RLB = new ArrayXd[bkgCount+1];

	ArrayXd* SNR = new ArrayXd[bkgCount+2];	
    
    //for average mean and cov
    
    MeansB[bkgCount].resize(SensorDims);
    CovsB[bkgCount].resize(SensorDims,SensorDims);  
    MeansB[bkgCount].setZero();
    CovsB[bkgCount].setZero();
    
    ArrayXd MeansT(SensorDims);
    ArrayXd MeansTinB(SensorDims);
    ArrayXd LT(SensorDims);
    MatrixXd CovsT(SensorDims,SensorDims);
    MatrixXd CovsLT(SensorDims,SensorDims);
    ArrayXd RLT(SensorDims);
    ArrayXd NoiseT(SensorDims);
    
    ArrayXd Llow(SensorDims);
    ArrayXd Lhigh(SensorDims);
    
    MatrixXd WlgthChosen(FeatureDims,1);
    MatrixXd BWChosen(FeatureDims,1);
    
    WlgthChosen = FeaturizeArray(Sen.getWlgth(),FeatureMtrx);
    BWChosen = FeaturizeArray(Sen.getBW(),FeatureMtrx);
    
//    WlgthChosen = FeatureMtrx*Sen.Wlgth.matrix();
    BWChosen = FeatureMtrx*Sen.getBW().matrix();

    ArrayXd WavelengthsSourceTemp, V_MeansSourceTemp;
    MatrixXd M_CovsSourceTemp;
  
//background 1,2,3, avg
	////B2001_15.BIN predefined 
    string pathtape7(m_FileDir->getmodtran()+"/tape7");
	string pathtape5(m_FileDir->getmodtran()+"/tape5");
    genTape5(Sen, Sce, -1, "DATA/B2001_15.BIN",pathtape5.c_str());  

	for (int i = 0; i<bkgCount;i++)
	{
		MeansB[i].resize(SensorDims);
		CovsB[i].resize(SensorDims,SensorDims); 
		string temppath = m_FileDir->getreflectance() + "/" + Bkg.getbackname()[i];       
		readWavelengthsMeansandCov(temppath.c_str(), WavelengthsSourceTemp, V_MeansSourceTemp, M_CovsSourceTemp);   
		LinearInterpolation1D(WavelengthsSourceTemp, V_MeansSourceTemp, Sen.getWlgth(), MeansB[i]); 
		LinearInterpolation2D(WavelengthsSourceTemp, M_CovsSourceTemp, Sen.getWlgth(), CovsB[i]);
         
		MeansB[bkgCount] = MeansB[bkgCount] + (MeansB[i]*Bkg.getbackperc()[i]/100);
		CovsB[bkgCount] += (CovsB[i]+(MeansB[i].matrix())*(MeansB[i].matrix().transpose()))*Bkg.getbackperc()[i]/100;        

		logfile(Writespec_albFile(Sen.getWlgth(),MeansB[i],m_FileDir),"//Writespec_albFile background"); 

		////linux may diff
		SetCurrentDirectoryA(m_FileDir->getmodtran().c_str());
		system((m_FileDir->getmodtran() + "/Mod4V3R1.exe").c_str());
//        cout<<"MeansB["<<i<<"]"<<" is running..."<<endl;
		logfile(i,"//MeansB[i] done");
		logfile(i,"ReadAllColumnFromTAPE7");
		ReadAllColumnFromTAPE7(pathtape7.c_str(),WVNo[i],SurfEmis[i],SolScat[i],GndRflt[i]);
		logfile(i,"ReadAllColumnFromTAPE7 done");
	}

    CovsB[bkgCount] = CovsB[bkgCount] - MeansB[bkgCount].matrix()*MeansB[bkgCount].matrix().transpose();
    logfile(Writespec_albFile(Sen.getWlgth(),MeansB[bkgCount],m_FileDir),"//Writespec_albFile background mean"); 
	if(bkgCount>1)
	{
		system((m_FileDir->getmodtran() + "/Mod4V3R1.exe").c_str());
	//    cout<<"MeansBavg"<<" is running..."<<endl; 
		logfile("avg","//MeansBavg done");
		ReadAllColumnFromTAPE7(pathtape7.c_str(),WVNo[bkgCount],SurfEmis[bkgCount],SolScat[bkgCount],GndRflt[bkgCount]);
	}
	else
	{
		WVNo[bkgCount] = WVNo[0];
		SurfEmis[bkgCount] = SurfEmis[0];
		SolScat[bkgCount] = SolScat[0];
		GndRflt[bkgCount] = GndRflt[0];
	}
//target 
    string temppath = m_FileDir->getreflectance() + "/" + Tgt.gettargetname();  
    readWavelengthsMeansandCov(temppath.c_str(), WavelengthsSourceTemp, V_MeansSourceTemp, M_CovsSourceTemp);
    LinearInterpolation1D(WavelengthsSourceTemp, V_MeansSourceTemp, Sen.getWlgth(), MeansT); 
    LinearInterpolation2D(WavelengthsSourceTemp, M_CovsSourceTemp, Sen.getWlgth(), CovsT);
    
    MeansTinB = MeansT * Tgt.gettargperc()/100 +(1-Tgt.gettargperc()/100)* MeansB[Tgt.gettarginback()-1];

	logfile(Tgt.gettarginback(),"//target in black NO."); 

    logfile(Writespec_albFile(Sen.getWlgth(),MeansTinB,m_FileDir),"//Writespec_albFile T in B");  
	////call modtran, different in LINUX
    system((m_FileDir->getmodtran() + "/Mod4V3R1.exe").c_str());
	logfile(bkgCount+1,"ReadAllColumnFromTAPE7");
    ReadAllColumnFromTAPE7(pathtape7.c_str(),WVNo[bkgCount+1],SurfEmis[bkgCount+1],SolScat[bkgCount+1],GndRflt[bkgCount+1]); 
	logfile(bkgCount+1,"ReadAllColumnFromTAPE7 done");

//fLBT<<MeansT<<endl;    
// for albedo  = 1, 0, 0.1, 0.6//   double albedo[4]={1,0,0.1,0.6};

	vector<double> albedo(2);
	albedo[0] = 1;
	albedo[1] = 0;

    for(size_t i = 0; i<albedo.size(); i++)
    {
         genTape5(Sen, Sce, albedo[i], "DATA/B2001_15.BIN",pathtape5.c_str()); 
		 ////call modtran, different in LINUX
         system((m_FileDir->getmodtran() + "/Mod4V3R1.exe").c_str()); 
		 logfile(bkgCount+2+i,"ReadAllColumnFromTAPE7");
         ReadAllColumnFromTAPE7(pathtape7.c_str(), WVNo[bkgCount+2+i],SurfEmis[bkgCount+2+i],SolScat[bkgCount+2+i],GndRflt[bkgCount+2+i]); 
		 logfile(bkgCount+2+i,"ReadAllColumnFromTAPE7 done");
//         cout<<"Albedo = "<<albedo[i]<<" is running..."<<endl;
		 logfile(i,"Albedo[i] done");
    }

//gaussian waveform
logfile("GenRespFunction");
    for(int i=0; i < bkgCount+extrarun;i++) 
	        GenRespFunction(SurfEmis[i],SolScat[i],GndRflt[i],
		                   SurfEmisO[i],SolScatO[i],GndRfltO[i],
		                   WVNo[0], Sen.getBW(), Sen.getWlgth()/1000);
logfile("GenRespFunction done");
    delete []WVNo;
//	delete []Trans;
	delete []SurfEmis;
	delete []SolScat;
	delete []GndRflt;
//	delete []DrctRflt;
//	delete []TotRad;

    MatrixXd CovBavePath =  WriteDiagonalMatrix(SolScatO[bkgCount+2] - SolScatO[bkgCount+3]) * 
                            Bkg.getbackscale() * (CovsB[bkgCount]/10000) * 
                            WriteDiagonalMatrix(SolScatO[bkgCount+2] - SolScatO[bkgCount+3]);                      

    MatrixXd Ddown(SensorDims,SensorDims);
    MatrixXd Tdown(SensorDims,SensorDims);
    
//    Llow = GndRfltO[bkgCount+4]+SolScatO[bkgCount];
//    Lhigh = GndRfltO[bkgCount+5]+SolScatO[bkgCount];


	Llow =  0.1*GndRfltO[bkgCount+2]+SolScatO[bkgCount];
	Lhigh = 0.6*GndRfltO[bkgCount+2]+SolScatO[bkgCount]; 	    

    ArrayXd ffactor = fixedfactor(Sen.getifov(), Sen.getaperture(), Sen.getOPTtrans(), Sen.getQeff(), Sen.gettint());
    ArrayXd fsigma_nbe_s = Sigma_nbe_s(Sen.getnbits(),Sen.getLSB(),Sen.getber());
    ArrayXd fsigma_nq_s = Sigma_nq_s(Sen.getLSB());
    
//cout<<"CovBavePath Done"<<endl;
logfile("CovBavePath Done");


    for(int i=0; i < bkgCount+1;i++)
    {
        CovsLB[i].resize(SensorDims,SensorDims);

        for(int j=0;j<SensorDims;j++)
              if (MeansB[i][j]<0.1) GndRfltO[i][j]=0;         
        
        Ddown = WriteDiagonalMatrix(GndRfltO[i]/(MeansB[i]/100));

        CovsLB[i] =  Bkg.getbackscale() * Ddown * CovsB[i]/10000 *  Ddown + CovBavePath;   

		LB[i] = GndRfltO[i] + SurfEmisO[i] + SolScatO[bkgCount];
		//Gndflt = downwelled + direct reflected in each case
		//Surfemis = Emission as a function of its temperature and emissivity, almost all zero in this case
		//SolScatO[bkgCount] = path radiance(upwelled radiance) in case "background average"
        RLB[i].resize(SensorDims);
        RLB[i]= ELM(LB[i],Llow,Lhigh);
        SNR[i].resize(SensorDims);
        
        NoiseB[i] = SensorProcessNoise(SNR[i],Sen.getFNoise(), ffactor,Sen.getrelcal(),Sen.getgainfac(),
              fsigma_nbe_s,fsigma_nq_s,LB[i],Sen.getWlgth()/1000,Sen.getBW(), Sen.getDegrade()); 
            
        CovsLB[i] = CovsLB[i] + WriteDiagonalMatrix( NoiseB[i]);    
        fLBT<<LB[i]<<endl; 
        fSNR<<SNR[i]<<endl; 
              
    }

    cout<<"Adding Noise Done"<<endl;
    logfile("Adding Noise Done");
 
	Tdown = WriteDiagonalMatrix(GndRfltO[Tgt.gettarginback()-1]/(MeansB[Tgt.gettarginback()-1]/100)); 

    CovsLT = (Tgt.gettargperc()/100*Tgt.gettargperc()/100)*Tgt.gettargscale() * Tdown * CovsT/10000 * Tdown +
             (1-Tgt.gettargperc()/100)*(1-Tgt.gettargperc()/100) * (CovsLB[Tgt.gettarginback()-1]-CovBavePath) + CovBavePath; 

    LT =  GndRfltO[bkgCount+1]+ SurfEmisO[bkgCount+1] + SolScatO[bkgCount]; 
    RLT = ELM(LT,Llow,Lhigh);
  
    SNR[bkgCount+1].resize(SensorDims);
    NoiseT = SensorProcessNoise(SNR[bkgCount+1],Sen.getFNoise(), ffactor,Sen.getrelcal(),Sen.getgainfac(),
              fsigma_nbe_s,fsigma_nq_s,LT,Sen.getWlgth()/1000,Sen.getBW(), Sen.getDegrade()); 
	CovsLT = CovsLT + WriteDiagonalMatrix(NoiseT);
	fLBT<<LT<<endl; 
    fSNR<<SNR[bkgCount+1]<<endl;     


    cout<<"LT, SNR output Done"<<endl;
	logfile("LT, SNR output Done");
    ArrayXd temp = ArrayXd::Ones(SensorDims)*0.5;
    MatrixXd coeffMatrix = (temp/(Lhigh-Llow)).matrix()*(temp/(Lhigh-Llow)).matrix().transpose();

    MatrixXd CovsRBaveInv = FeaturizeMatrix(((CovsLB[bkgCount].array() * coeffMatrix.array()).matrix()),FeatureMtrx).inverse();

    MatrixXd low = FeaturizeArray((MeansT/100-RLB[bkgCount]),FeatureMtrx).matrix().transpose() * 
                    CovsRBaveInv * 
                   FeaturizeArray((MeansT/100-RLB[bkgCount]),FeatureMtrx).matrix();
    MatrixXd w = CovsRBaveInv * (FeaturizeArray((MeansT/100-RLB[bkgCount]),FeatureMtrx).matrix())/low(0);
     
//w<<-0.261008420000000,-0.0176464770000000,-0.0326333900000000,0.225798510000000,0.472442750000000,0.817609030000000,0.597607360000000,0.154493870000000,-0.418021830000000,-0.349629320000000,-0.378337180000000,-0.392110630000000,-0.0824135840000000,0.0567781340000000,0.439369290000000,-0.616467900000000,-0.926037180000000,0.432299440000000,0.0804048880000000,-0.145172300000000,-0.109771410000000,-0.0366633530000000,-0.0232091520000000,-0.00985189440000000,0.0105217040000000,0.00729003660000000,0.0365207610000000,-0.000620496800000000,-0.0264728250000000,-0.0439757160000000,-0.0421326610000000,0.0132648150000000,0.0653972470000000,0.0497869050000000,-0.0270363870000000,-0.00982702770000000,0.0437304540000000,-0.0348783000000000,0.00325514730000000,-0.0214939340000000,-0.0462768050000000,-0.0272654370000000,0.0153502060000000,0.0545333650000000,0.113224690000000,0.148193630000000,-0.0121015070000000,0.0851857210000000,0.0283402140000000,0.154458440000000,0.0843239120000000,0.0878276760000000,0.0513694380000000,-0.0336532560000000,-0.124117050000000,-0.143684060000000,-0.0720325750000000,-0.0518801460000000,-0.0591967910000000,-0.0560700840000000,0.0927735370000000,-0.190341810000000,1.06215240000000e-12,4.26069840000000e-06,0.00138633540000000,-0.00277829560000000,0.0281898230000000,0.226105440000000,0.0155218250000000,-0.0239099360000000,0.222584890000000,0.198467470000000,0.0305272880000000,-0.102863860000000,-0.318593130000000,0.211679280000000,-0.00625364600000000,-0.180031790000000,0.0214184100000000,-0.167251820000000,-0.0572059910000000,-0.0650003880000000,-0.0982796440000000,0.0220583970000000,-0.0172472160000000,0.0285595370000000,0.0150018580000000,-0.00666273400000000,0.0364149130000000,0.111319960000000,0.0796815490000000,0.00847889750000000,-0.00443656240000000,-0.0288518700000000,-0.0632022180000000,-0.197982520000000,1.15453290000000,-0.278354410000000,-0.304048200000000,-0.163407070000000,0.337927610000000,0.399878130000000,0.705311170000000,-0.571593300000000,-0.382312640000000,-0.431613770000000,-0.0948073150000000,0.687234650000000,1.72519210000000,1.27941400000000,-0.910963130000000,-1.41471180000000,-0.698180300000000,0.0233450000000000,0.212713280000000,0.131649750000000,0.0591320450000000,-0.0644962440000000,0.0557122100000000,-0.0808078290000000,-0.00169258060000000,0.122156630000000,0.0898191270000000,0.0429405680000000,-0.00220469670000000,-0.000116832160000000;    


//cout<<"w calculation Done"<<endl;
	logfile(low(0),"low");
	logfile("w calculation Done");
#ifdef DEBUG
//fdebug<<endl<<"CovsLB[1]"<<endl<<CovsLB[1]<<endl<<endl;
//fdebug<<endl<<"CovsLB[0]"<<endl<<CovsLB[0]<<endl<<endl;
//fdebug<<endl<<"CovsLT"<<endl<<CovsLT<<endl<<endl;


//fdebug<<endl<<"MeansTinB"<<endl<<MeansTinB<<endl<<endl;
//fdebug<<endl<<"MeansB[0]"<<endl<<MeansB[0]<<endl<<endl;
//fdebug<<endl<<"MeansT"<<endl<<MeansT<<endl<<endl;

//fdebug<<endl<<"coeffMatrix"<<endl<<coeffMatrix<<endl<<endl;  
//fdebug<<endl<<"Llow"<<endl<<Llow<<endl<<endl;
//fdebug<<endl<<"Lhigh"<<endl<<Lhigh<<endl<<endl;
//fdebug<<endl<<"LB1"<<endl<<LB[0]<<endl<<endl;

//fdebug<<endl<<"CovsB[1]"<<endl<<CovsB[1]<<endl<<endl;
//fdebug<<endl<<"MeansBgrass"<<endl<<MeansB[1]<<endl<<endl;
//fdebug<<endl<<"CovsLB[1]"<<endl<<CovsLB[1]<<endl<<endl;

//fdebug<<endl<<"CovBavePath"<<endl<<CovBavePath<<endl<<endl;
//fdebug<<endl<<"CovsT/10000"<<endl<<CovsT/10000 <<endl<<endl;  
//fdebug<<endl<<"CovsTS"<<endl<<WriteDiagonalMatrix(TransO[bkgCount+1])*Tdown * CovsT/10000 * Tdown*WriteDiagonalMatrix(TransO[bkgCount+1])<<endl<<endl;
//fdebug<<endl<<"CovsRBave"<<endl<<(CovsLB[bkgCount].array() * coeffMatrix.array())<<endl<<endl;
//fdebug<<endl<<"bscov_targinback"<<endl<<(CovsLB[Tgt.gettarginback()-1]-CovBavePath)<<endl<<endl;
//fdebug<<endl<<"CovsRB0"<<endl<<(CovsLB[0].array() * coeffMatrix.array())<<endl<<endl;
//fdebug<<endl<<"CovsRB1"<<endl<<(CovsLB[1].array() * coeffMatrix.array())<<endl<<endl;
//fdebug<<endl<<"CovsRB2"<<endl<<(CovsLB[2].array() * coeffMatrix.array())<<endl<<endl;
//fdebug<<endl<<"CovsRT"<<endl<<(CovsLT.array() * coeffMatrix.array())<<endl<<endl;
//fdebug<<endl<<"CovsLT"<<endl<<CovsLT.array()<<endl<<endl;
//fdebug<<endl<<"CovsRBaveInv"<<endl<<CovsRBaveInv<<endl<<endl;
//fdebug<<endl<<"coeffMatrix"<<endl<<coeffMatrix<<endl<<endl;    fdebug.close();

//fdebug<<endl<<"RLB[bkgCount]"<<endl<<RLB[bkgCount]<<endl<<endl;
//fdebug<<endl<<"RLB[0]"<<endl<<RLB[0]<<endl<<endl;
//fdebug<<endl<<"RLT"<<endl<<RLT<<endl<<endl;
//fdebug<<endl<<"WlgthChosen"<<endl<<WlgthChosen<<endl<<endl;

//fdebug<<endl<<"low(0)"<<endl<<low(0)<<endl<<endl; 
//fdebug<<endl<<"w"<<endl<<w<<endl<<endl;  
#endif

//	delete []TransO;
	delete []SurfEmisO;
	delete []SolScatO;
	delete []GndRfltO;
//	delete []DrctRfltO;
//	delete []TotRadO;

    MatrixXd mt = w.transpose()*(FeatureMtrx*((RLT-RLB[bkgCount]).matrix()));
    cout<<"MT="<<mt<<endl;
	logfile(mt(0),"MT");
    MatrixXd ct = w.transpose()*FeaturizeMatrix((CovsLT.array()*coeffMatrix.array()).matrix(),FeatureMtrx)*w;  

    cout<<"ct="<<sqrt(ct(0))<<endl;
    logfile(ct(0),"CT");
    MatrixXd* mb = new MatrixXd[bkgCount];
    MatrixXd* cb = new MatrixXd[bkgCount];
    
    ArrayXd meansROC(bkgCount+1);
    ArrayXd stdROC(bkgCount+1);
    for(int i=0;i<bkgCount;i++)
    {
        mb[i] = w.transpose()*(FeatureMtrx*(RLB[i]-RLB[bkgCount]).matrix());
        cb[i] = w.transpose()*FeaturizeMatrix((CovsLB[i].array()*coeffMatrix.array()).matrix(),FeatureMtrx)*w;
        cout<<"mb["<<i<<"]="<<mb[i]<<endl;
        cout<<"cb["<<i<<"]="<<sqrt(cb[i](0))<<endl;
		logfile(mb[i],"MBi");
		logfile(cb[i],"CBi");
        meansROC[i]=mb[i](0);
        stdROC[i]=sqrt(cb[i](0));
    }

    meansROC[bkgCount]=mt(0);
    stdROC[bkgCount]=sqrt(ct(0));

//    mm[0]=0;
//    dd[0]=0.0093;
//    mm[bkgCount]=0.0488662;
//    dd[bkgCount]=0.0104;

	delete []MeansB;
    delete []CovsB;
    delete []CovsLB;
	delete []NoiseB;
	delete []LB;
	delete []RLB;

	delete []mb;
	delete []cb;
    
    cout<<meansROC<<endl;
    cout<<stdROC<<endl;


    ArrayXd outPdPfa=ROCgen(meansROC, stdROC, Bkg.getbackperc(), Pro.getpfa(),bkgCount);

	fROC<<outPdPfa;

    int p = 0;
    while(p!=Pro.getpfa().size()) cout<<Pro.getpfa()[p++]<<endl;


    fLBT.close();
    fSNR.close();
    fROC.close();
    fdebug.close();

	flog<<endl<<endl<<endl;
	flog.close();
    

    delete m_FileDir;
//    system("PAUSE");
  
    return EXIT_SUCCESS;
	}
	catch(GeneralException)
	{
		return -1;
	}
}

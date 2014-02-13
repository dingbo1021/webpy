#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <vector>
#include <string>
#include <iomanip>  
#include <Eigen/Dense>
#include <stdio.h>
#include <windows.h>
#include <math.h>
#include <errno.h>

#include <time.h> 

using namespace std;
using namespace Eigen;


inline double ELM(double L, double L1, double L2) {return 0.1f+0.5f*(L-L1)/(L2-L1);};
inline ArrayXd ELM(VectorXd L, VectorXd L1, VectorXd L2) {return 0.1f+0.5f*(L.array()-L1.array())/(L2.array()-L1.array());};
inline ArrayXd ELM(ArrayXd L, ArrayXd L1, ArrayXd L2) {return 0.1+0.5*(L-L1)/(L2-L1);};
inline ArrayXd ELM_inverse(VectorXd r, VectorXd L1, VectorXd L2) {return 2*r.array()*(L2.array()-L1.array())-0.2f*L2.array()+1.2f*L1.array();};
inline ArrayXd FeaturizeArray(ArrayXd X, MatrixXd FeatureMtrx){return (FeatureMtrx*X.matrix()).array();}
inline MatrixXd FeaturizeMatrix(MatrixXd X, MatrixXd FeatureMtrx){return (FeatureMtrx*X*FeatureMtrx.transpose());}

class GeneralException {};
////log file might get too large
ofstream flog("C:/webpy/systemlog.log",ios::app);

template <class T>
void logfile(T input, char* comment = " ")
{
	static int logCount = 0;
	if(0==logCount)
	{
		logCount++;
		if(!flog.is_open())
		{
			logfile("Logfile cannot be established.");
			throw GeneralException();
		}
	}
	time_t t = time(0); 
	char tmp[64]; 
	strftime( tmp, sizeof(tmp), "%Y/%m/%d %X",localtime(&t) ); 
	flog<<tmp<<"   "<<input<<"    "<<comment<<endl;
}

class FileDir
{
	public:
		FileDir();
		FileDir(const char* inputpathconfig);
		FileDir(char* pfile, char* sensor, char* reflectance,  char* modtranin);    

		string getpfile() const{return pfile;};
		string getsensor() const{return sensor;};
		string getreflectance() const{return reflectance;};
		string getmodtran() const{return modtran;};

		const char* getLBT() const{return LBT.c_str();};
		const char* getSNR() const{return SNR.c_str();};
		const char* getROC() const{return ROC.c_str();};


		const char* getForDebug() const{return ForDebug.c_str();};
		const char* getSystemLog() const{return SystemLog.c_str();};

		const char* getwavelengthchosen() const{return wavelengthchosen.c_str();};


		void showall() const;
		void logall() const;
       
	private:
		const char* pathconfig;
		string pfile;
		string sensor;
		string reflectance;
		string modtran;
		string LBT;
		string SNR;
		string ROC;

		string ForDebug;
		string SystemLog;

		string wavelengthchosen;
  
};
////default path for those files
FileDir::FileDir(const char* inputpathconfig):pathconfig(inputpathconfig)
{
    ifstream fin(inputpathconfig);
	if(!fin.is_open())
	{
		logfile("FileDir setting wrong, check pathconfig.ini.");
		throw GeneralException();
	}
	getline(fin,pfile);
	getline(fin,sensor);
	getline(fin,reflectance);
	getline(fin,modtran);

	getline(fin,LBT);
	getline(fin,SNR);
	getline(fin,ROC);

	getline(fin,ForDebug);
	getline(fin,SystemLog);

	getline(fin,wavelengthchosen);

	fin.close(); 
	logfile("FileDir setting completed.");
}
////predefined pathconfig file path
FileDir::FileDir():pfile("C:/FASSP_EMMETT/v2.4/pfiles"),sensor("C:/FASSP_EMMETT/v2.4/sens"), 
	reflectance("C:/FASSP_EMMETT/v2.4/refl"),modtran("D:/RIT/JohnThesis/MOD4v3r1")
{
    ;
}
//set path for those files
FileDir::FileDir(char* pfilein, char* sensorin, char* reflectancein, char* modtranin)
{
    if(pfilein&&sensorin&&reflectancein)  
    {
        pfile = pfilein;                            
        sensor = sensorin;                            
        reflectance = reflectancein; 
        modtran = modtranin;                            
    }
	else
		throw GeneralException();

}

void FileDir::showall() const
{
     cout<<pfile<<endl;
     cout<<sensor<<endl;
     cout<<reflectance<<endl;
     cout<<modtran<<endl;   
}

void FileDir::logall() const
{
	logfile(pfile); 
	logfile(sensor); 
	logfile(reflectance); 
	logfile(modtran); 

}


bool Writespec_albFile(ArrayXd& Wavelengths, ArrayXd& Reflectances, FileDir* m_FileDir)
{
	if(Wavelengths.size()!=Reflectances.size()) 
	{
		logfile("Writespec_albFile Wavelengths.size()!=Reflectances.size()","/ERROR");
		throw GeneralException();
	}
	//// // or \\ in linux
	string albpath = m_FileDir->getmodtran() + "/DATA/spec_alb.dat";
    ofstream fouttemp(albpath.c_str());
	if(!fouttemp.is_open())
	{
		logfile("Writespec_albFile cannot write","/ERROR");
		throw GeneralException();	
	}
	fouttemp<<"1 alb file "<<endl;
	for(int i=0; i< Wavelengths.size(); i++)
		fouttemp<<"   "<<Wavelengths(i)/1000<<" "<<Reflectances(i)/100<<endl;
	fouttemp<<"!"<<endl;
	fouttemp.close();
	return true;

}

void ReadAllColumnFromTAPE7(const char* tape7path, ArrayXd& WVNo, ArrayXd& SurfEmis, ArrayXd& SolScat, ArrayXd& GndRflt)
{
	int i=0;
	VectorXd dataOut;
	logfile(tape7path);
	static int readCount = 1;
	if(!tape7path)
	{
		logfile(readCount,"//ReadAllColumnFromTAPE7 tape7path NULL");
		throw GeneralException();
	}
    ifstream fin(tape7path);
    string s;
	if(!fin.is_open())	
	{
		logfile(readCount,"//ReadAllColumnFromTAPE7 tape7path FAILS");
		throw GeneralException();
	}
	streamoff keyposition=0,offset=0;
	for(int j=0;j<8;j++) getline(fin,s); //targets are in line 9 
	fin>>s;
	double Llmts = (double)atof(s.c_str());
	fin>>s;
	double Ulmts = (double)atof(s.c_str());
	fin>>s;
	double Step = (double)atof(s.c_str());

	int size = (int)((Ulmts-Llmts)/Step+1.999);
	WVNo.resize(size);
//	Trans.resize(size);
	SurfEmis.resize(size);
	SolScat.resize(size);
	GndRflt.resize(size);
//	DrctRflt.resize(size);
//	TotRad.resize(size);
	float temp[7];

	for(int j=0;j<3;j++) getline(fin,s);

	while(i<size)
	{
		getline(fin,s);
		////bo  sscanf_s ???
		sscanf(s.c_str(),
			"%8f %11f %*11f %*11f %11f %11f %*11f %11f %11f %f %*f %*f %*f",
			&temp[0],&temp[1],&temp[2],&temp[3],&temp[4],&temp[5],&temp[6]);

		WVNo(i) = temp[0];
//		Trans(i) = temp[1];
		SurfEmis(i) = temp[2];
		SolScat(i) = temp[3];
		GndRflt(i) = temp[4];
//		DrctRflt(i) = temp[5];
//		TotRad(i) = temp[6];
		i++;
	}  
	fin.close();
	logfile(readCount,"//ReadAllColumnFromTAPE7");

	readCount++;


}

void GenRespFunction(ArrayXd dataIn1,  ArrayXd dataIn2,  ArrayXd dataIn3,  
	                 ArrayXd &dataOut1,ArrayXd &dataOut2,ArrayXd &dataOut3,
				     ArrayXd dataWVNo, ArrayXd bandwidth, ArrayXd SensorWlength)
{
	static int genRespCount = 1;
	int WNoCnt  = dataWVNo.size();
	int WNoStep = (int)(dataWVNo(1) - dataWVNo(0));
	int sensorsize = SensorWlength.size();
	
	double UlimitWN = dataWVNo(WNoCnt-1);
	double LlimitWN = dataWVNo(0);
	double SQRT2ln2 = 1.17741002251547469;
	double SQRT2Pi = 2.506628274631;
	int index1 = 0, index2 = 0, Seglength = 0;
	
	ArrayXd AllOne = ArrayXd::Ones(SensorWlength.size());
	ArrayXd dataFreq = ArrayXd::Ones(WNoCnt)*10000/dataWVNo;
	ArrayXd WaveNoTOP = AllOne*10000/(SensorWlength - 3*bandwidth) ;
	ArrayXd WaveNoButm = AllOne*10000/(SensorWlength + 3*bandwidth);
	ArrayXd Sigma = bandwidth*0.5/SQRT2ln2;
	ArrayXd factor =  AllOne/(Sigma*SQRT2Pi);
	ArrayXd DataSeg1,DataSeg2,DataSeg3;

	ArrayXd WlengthSeg;
	ArrayXd WVNoSeg;
	ArrayXd SigmaSeg;
	ArrayXd Temp;
	
	dataOut1.resize(sensorsize);
	dataOut2.resize(sensorsize);
	dataOut3.resize(sensorsize);
	//dataOut4.resize(sensorsize);
	//dataOut5.resize(sensorsize);
	//dataOut6.resize(sensorsize);

	for(int i=0; i < SensorWlength.size();i++)
	{
		index1 = (int)(WaveNoButm(i)-LlimitWN)/WNoStep+1;
		index2 = (int)(WaveNoTOP(i)-LlimitWN)/WNoStep;
		Seglength = index2 - index1 + 1;
//cout<<index1<<"  "<<index2<<endl;
//cout<<LlimitWN<<"  "<<UlimitWN<<" "<<WNoCnt<<endl;
//system("PAUSE");
		DataSeg1.resize(Seglength);
		DataSeg2.resize(Seglength);
		DataSeg3.resize(Seglength);
		//DataSeg4.resize(Seglength);
		//DataSeg5.resize(Seglength);
		//DataSeg6.resize(Seglength);
		WVNoSeg = dataWVNo.segment(index1,Seglength);
		DataSeg1 = dataIn1.segment(index1,Seglength);
		DataSeg2 = dataIn2.segment(index1,Seglength);
		DataSeg3 = dataIn3.segment(index1,Seglength);
		//DataSeg4 = dataIn4.segment(index1,Seglength);
		//DataSeg5 = dataIn5.segment(index1,Seglength);
		//DataSeg6 = dataIn6.segment(index1,Seglength);

		WlengthSeg = 1e4*(ArrayXd::Ones(Seglength))/WVNoSeg;

		DataSeg1 = 0.1*DataSeg1*WVNoSeg*WVNoSeg;
		DataSeg2 = 0.1*DataSeg2*WVNoSeg*WVNoSeg;
		DataSeg3 = 0.1*DataSeg3*WVNoSeg*WVNoSeg;
		//DataSeg4 = 0.1*DataSeg4*WVNoSeg*WVNoSeg;
		//DataSeg5 = 0.1*DataSeg5*WVNoSeg*WVNoSeg;

		Temp = factor(i)*exp(-(WlengthSeg-SensorWlength(i))*(WlengthSeg-SensorWlength(i))/(2*Sigma(i)*Sigma(i)));
		Temp = Temp/Temp.sum();

		dataOut1(i)=(Temp*DataSeg1).sum();
		dataOut2(i)=(Temp*DataSeg2).sum();
		dataOut3(i)=(Temp*DataSeg3).sum();
		//dataOut4(i)=(Temp*DataSeg4).sum();
		//dataOut5(i)=(Temp*DataSeg5).sum();
		//dataOut6(i)=(Temp*DataSeg6).sum();

	}
	logfile(genRespCount,"//genResp done");
	genRespCount++;
}

void parsemultiinput(string s, vector<double> &output)
{
    string subs;
    const char *c = s.c_str();
	int index = 0;
	for(size_t i = 0; i < strlen(c); i++)
	{
        if(*(c + i) == ','||*(c + i) == '/')
        {
            if(index == 0) 
            {
                   index = i+1;
                   continue; 
            }
            subs=s.substr(index,i-index);     
//            cout<<subs<<"~double"<<endl;                               
            index = i+1;
            output.push_back(atof(subs.c_str()));
        }
    }
	if(output.empty())
	{
		logfile("parsemultiinput vector<double> is EMPTY","/ERROR");
		throw GeneralException();
	}
}
void parsemultiinput(string s, vector<int> &output)
{
    string subs;
    const char *c = s.c_str();
	int index = 0;
	for(size_t i = 0; i < strlen(c); i++)
	{
        if(*(c + i) == ','||*(c + i) == '/')
        {
            if(index == 0) 
            {
                   index = i+1;
                   continue; 
            }
            subs=s.substr(index,i-index);     
//            cout<<subs<<"~double"<<endl;                               
            index = i+1;
            output.push_back(atoi(subs.c_str()));
        }
    }
	if(output.empty())
	{
		logfile("parsemultiinput vector<int> is EMPTY","/ERROR");
		throw GeneralException();
	}
}
void parsemultiinput(string s, vector<string> &output)
{
    string subs;   
    const char *c = s.c_str();
    int index = 0;
    for(size_t i = 0; i < strlen(c); i++)
    {
    
        if(*(c + i) == ','||*(c + i) == '/')
        {
            if(index == 0) 
            {
                   index = i+1;
                   continue; 
            }
            subs=s.substr(index,i-index);     
//            cout<<subs<<"~string"<<endl;                               
            index = i+1;
            output.push_back(subs);
        }
    }
	if(output.empty())
	{
		logfile("parsemultiinput vector<string> is EMPTY","/ERROR");
		throw GeneralException();
	}
}

void parsesingle(string s, char* buffer)
{
	if(buffer)
	{
		if(s.c_str()[0]=='/')
			sscanf(s.c_str(),"%*[/] %[0-9.] %*[/]",buffer);
		else
			strcpy(buffer, s.c_str());
	}
	else
	{
			
		logfile("parsesingle buffer is NULL","/ERROR");
		throw GeneralException();
	}
}

void cutcomment(string &s)
{
    string subs;
    const char *c = s.c_str();
	int index = 0;
	for(size_t i = 0; i < strlen(c); i++)
	{
        if(*(c + i) == ';')
        {
            subs=s.substr(index,i-index);     
            s = subs + ",";  
//            cout<<s<<"~cutcomment"<<endl;                       
        }
    } 
}

class Target
{
public:
	Target(const char* pfilepath);
	string gettargetname() const{return targetname;};
	double gettargscale() const{return targscale;};
	int gettnsamp() const{return tnsamp;};
	double gettargperc() const{return targperc;};
	int gettarginback() const{return targinback;};
    
private:
	string targetname;
	double targscale;
	int tnsamp;
	double targperc;
	int targinback; 

};
Target::Target(const char* pfilepath)
{
	if(!pfilepath)
	{
		logfile("Target construct fails, pfilepath NULL","/ERROR");
		throw GeneralException();
	}
	ifstream fin(pfilepath);
	if(!fin.is_open())
	{
		logfile("Target construct fails, pfilepath wrong","/ERROR");
		throw GeneralException();
	}
    string s;
	char buffer[100];
	while(fin >> s)
	{
		if(s=="targetname")
		{	
            fin>>s;
			targetname = s;
			getline(fin,s);
			fin>>s;
			fin>>s;
			targscale = (double)atof(s.c_str());
			getline(fin,s);
			fin>>s;
			fin>>s;
			tnsamp = atoi(s.c_str());
			getline(fin,s);
			fin>>s;
			fin>>s;
			parsesingle(s,buffer);
			targperc = atof(buffer);
			getline(fin,s);
			fin>>s;
			fin>>s;
			targinback = atoi(s.c_str());
        }
	}
	logfile("Target construct over");

}



class Modtran
{
    public:
        Modtran();   
        Modtran(FileDir* FileDir);
        int runcode;

    private:
		const char* dir;
        double v1;
        double v2;
        int res1;
        int res2;
            
};
////predefined path
Modtran::Modtran():dir("D:/RIT/JohnThesis/MOD4v3r1"),runcode(1),v1(0),v2(0),res1(15),res2(20)
{
;              
}

Modtran::Modtran(FileDir* FileDir):dir(FileDir->getmodtran().c_str()),runcode(1),v1(0),v2(0),res1(15),res2(20)
{
;                 
}

class Background
{
public:
    Background(const char* pfilepath);
       
	int getnumback() const {return numback;};
	std::vector<string> getbackname() const {return backname;};
	std::vector<double> getbackperc() const {return backperc;};
	double getbackscale() const {return backscale;};
	int getbnsamp() const {return bnsamp;};
	int getbackdof() const {return backdof;};
       
private:
	int numback;
	std::vector<string> backname;
	std::vector<double> backperc;
	double backscale;
	int bnsamp;
	int backdof; 


};

Background::Background(const char* pfilepath)
{
	if(!pfilepath)
	{
		logfile("Background construct fails, pfilepath NULL","/ERROR");
		throw GeneralException();
	}
	ifstream fin(pfilepath);
	if(!fin.is_open())
	{
		logfile("Background construct fails, pfilepath WRONG","/ERROR");
		throw GeneralException();
	}

	string s,subs;
	int j = 0;
	int index = 1;
	while(fin >> s)
	{
		if(s=="backname")
		{

			getline(fin,s);           
			parsemultiinput(s,backname);

			fin>>s;
			fin>>s;               
			backscale = atof(s.c_str());
                
 			getline(fin,s);
			fin>>s;
			fin>>s;               
			bnsamp = atoi(s.c_str());   
                
 			getline(fin,s);             
			fin>>s;
			getline(fin,s); 
			
			parsemultiinput(s,backperc);
				
			fin>>s;
			fin>>s; 
       
			backdof = atoi(s.c_str());  
			getline(fin,s);
                
			fin>>s;
			fin>>s;
			numback = atoi(s.c_str()); 
			logfile("Background construct over");        

		}

	}


}

class Scene
{
public:
    Scene(const char* pfilepath);
	double getmetrange(){return metrange;};
	double getsolangle(){return solangle;};
	double getgndalt(){return gndalt;};
	int	   getmodel(){return model;};
	double getshadeperc(){return shadeperc;};
	double getskyperc(){return skyperc;};
	double getwss(){return wss;};
	int    getihaze(){return ihaze;};
	int    geticld(){return icld;};
	int    getruncode(){return runcode;};

	void    seticld(int input){icld = input;};

    void showall(void);	      
       
private:
	double      metrange;
	double      solangle;
    double      gndalt;
    int         model;
    double      shadeperc;
    double      skyperc;
    double      wss;
    int         ihaze;
    int         icld;
    int         runcode;


};

void Scene::showall(void)
{
     cout<<"metrange = "<<metrange<<endl;
     cout<<"solangle = "<<solangle<<endl;
     cout<<"gndalt = "<<gndalt<<endl;
     cout<<"model = "<<model<<endl;
     cout<<"shadeperc = "<<shadeperc<<endl;
     cout<<"skyperc = "<<skyperc<<endl;
     cout<<"wss = "<<wss<<endl;
     cout<<"ihaze = "<<ihaze<<endl;
     cout<<"icld = "<<icld<<endl;
     cout<<"runcode = "<<runcode<<endl;                        
}

Scene::Scene(const char* pfilepath)
{
	if(!pfilepath)
	{
		logfile("Scene construct fails, pfilepath NULL","/ERROR");
		throw GeneralException();
	}
	ifstream fin(pfilepath);
    string s,subs;
	if(!fin.is_open())
	{
		logfile("Scene construct fails, pfilepath Wrong","/ERROR");
		throw GeneralException();
	}
	while(fin >> s)
	{
		if(s=="metrange")
		{
			fin>>s;
			metrange = atof(s.c_str());
                
			getline(fin,s);
			getline(fin,s);
			fin>>s;
			fin>>s;
			solangle = atof(s.c_str()); 
              
 			getline(fin,s);
			fin>>s;
			fin>>s;               
			gndalt = atof(s.c_str());   
                
 			getline(fin,s);             
			fin>>s;
			fin>>s; 
			model = atoi(s.c_str()); 
			    
 			getline(fin,s);
			getline(fin,s);
			fin>>s;
			fin>>s; 
			shadeperc = atof(s.c_str());  
                
			getline(fin,s);
			fin>>s;
			fin>>s;
			skyperc = atoi(s.c_str()); 
         
			getline(fin,s);  
			fin>>s;
			fin>>s;
			wss = atof(s.c_str());
                
			getline(fin,s);  
			fin>>s;
			fin>>s;
			ihaze = atoi(s.c_str());                
                
			getline(fin,s);
			getline(fin,s);  
			fin>>s;
			fin>>s;
			icld = atoi(s.c_str());  
                              
			getline(fin,s);  
			fin>>s;
			fin>>s;
			runcode = atoi(s.c_str());                  
           
			logfile("Scene construct over");      

		}
	}
}

class Sensor
{
public:
    Sensor(const char* pfilepath, FileDir* FileDir);
       
	string getsensorfile() const{return sensorfile;};
	double getviewangle() const{return viewangle;};
	double getnoisefac() const{return noisefac;};
    double getgainfac() const{return gainfac;};
    double getrelcal() const{return relcal;};
    int    getnbits() const{return nbits;};
    double getber() const{return ber;};
    double getplatalt() const{return platalt;};
    double gettint() const{return tint;};
    
	int	   getndim() const{return ndim;};
    double getaperture() const{return aperture;};
    double getifov() const{return ifov;};
    double getelfnfac() const{return elfnfac;};
	int    getlwir() const{return lwir;};
        
    ArrayXd getWlgth() const{return Wlgth;};
    ArrayXd getBW() const{return BW;};
    ArrayXd getOPTtrans() const{return OPTtrans;};
    ArrayXd getQeff() const{return Qeff;};
    ArrayXd getLSB() const{return LSB;};
    ArrayXd getFNoise() const{return FNoise;};
    ArrayXd getDegrade() const{return Degrade;};
      
	int getwNohigh() const{return wNohigh;};
	int getwNolow() const{return wNolow;};    
	std::vector<int> getbackchosen() const {return backchosen;};
    void showall(void);	      
       
private:
	string sensorfile;
	double viewangle;
    double noisefac;
    double gainfac;
    double relcal;
    int    nbits;    
    double ber;    
    double platalt;
    double tint;
    
    int ndim;
    double aperture;
    double ifov;
    double elfnfac;
    int lwir;
	std::vector<int> backchosen;

    ArrayXd Wlgth;
    ArrayXd BW;
    ArrayXd OPTtrans;
    ArrayXd Qeff;
    ArrayXd LSB;
    ArrayXd FNoise;
    ArrayXd Degrade;
      
    int wNohigh;
    int wNolow;


};

void Sensor::showall(void)
{
     cout<<"sensorfile = "<<sensorfile<<endl;
     cout<<"viewangle = "<<viewangle<<endl;
     cout<<"noisefac = "<<noisefac<<endl;
     cout<<"gainfac = "<<gainfac<<endl;
     cout<<"relcal = "<<relcal<<endl;
     cout<<"nbits = "<<nbits<<endl;
     cout<<"ber = "<<ber<<endl;
     cout<<"platalt = "<<platalt<<endl;
     cout<<"aperture = "<<aperture<<endl;  
     cout<<"ifov = "<<ifov<<endl;
     cout<<"tint = "<<tint<<endl;        
     cout<<"OPTtrans[0~3] = "<<OPTtrans[0]<<" "<<OPTtrans[1]<<" "<<OPTtrans[2]<<endl;
     cout<<"Qeff[0~3] = "<<Qeff[0]<<" "<<Qeff[1]<<" "<<Qeff[2]<<endl; 
     cout<<"BW[0~3] = "<<BW[0]<<" "<<BW[1]<<" "<<BW[2]<<endl; 
     cout<<"Wv[0~3] = "<<Wlgth[0]<<" "<<Wlgth[1]<<" "<<Wlgth[2]<<endl;
     cout<<"Degrade[0~3] = "<<Degrade[0]<<" "<<Degrade[1]<<" "<<Degrade[2]<<endl; 
	 cout<<"Backchosen[0~3] = "<<backchosen[0]<<" "<<backchosen[1]<<" "<<backchosen[2]<<endl; 
}
Sensor::Sensor(const char* pfilepath, FileDir* FileDir)
{
	if(!(pfilepath&&FileDir))
	{
		logfile("Sensor construct fails, pfilepath or FileDir NULL","/ERROR");
		throw GeneralException();
	}
                    
	ifstream fin(pfilepath);
    string s,subs;
    float temp[7];

	if(!fin.is_open())
    {
		logfile("Sensor construct fails, pfilepath Wrong","/ERROR");
		throw GeneralException();
	}
	while(fin >> s)
	{
		if(s=="sensorfile")
		{
			fin>>s;
			sensorfile = s;
                
			getline(fin,s);
			fin>>s;
			fin>>s;
			viewangle = atof(s.c_str()); 
              
 			getline(fin,s);
			fin>>s;
			fin>>s;               
			noisefac = atof(s.c_str());   
                
 			getline(fin,s);             
			fin>>s;
			fin>>s; 
			gainfac = atof(s.c_str()); 
			    
 			getline(fin,s);
			fin>>s;
			fin>>s; 
			relcal = atof(s.c_str());  
                
			getline(fin,s);
			fin>>s;
			fin>>s;
			nbits = atoi(s.c_str()); 
         
			getline(fin,s);  
			fin>>s;
			fin>>s;
			ber = atof(s.c_str());
                
			getline(fin,s);  
			fin>>s;
			fin>>s;
			platalt = atof(s.c_str());                
                
			getline(fin,s);
			fin>>s;
			fin>>s;
			tint = atof(s.c_str());

			getline(fin,s);
			fin>>s;
			fin>>s;
			parsemultiinput(s,backchosen);
	    }
    }


    string sensorfilepath = FileDir->getsensor() + "/" + sensorfile;  
	ifstream fin2(sensorfilepath.c_str()); 
	if(!fin2.is_open())
	{
		logfile("Sensor file path Wrong","/ERROR");
		throw GeneralException();
	}
	while(fin2 >> s)
	{
		if(s=="ndim")
		{
//               system("PAUSE");
			fin2>>s;
			ndim = atoi(s.c_str());
                  
			getline(fin2,s);  
			fin2>>s; 
			fin2>>s;                               
			aperture = atof(s.c_str());
                
			getline(fin2,s);
			fin2>>s;
			fin2>>s;                
			ifov = atof(s.c_str());
                
			getline(fin2,s);
			fin2>>s;
			fin2>>s;                                
			elfnfac = atof(s.c_str());
                
			getline(fin2,s);
			fin2>>s;
			fin2>>s;                
			lwir = atoi(s.c_str());
                
			Wlgth.resize(ndim);
			BW.resize(ndim);
			OPTtrans.resize(ndim);
			Qeff.resize(ndim);
			LSB.resize(ndim);
			FNoise.resize(ndim);
			Degrade.resize(ndim);
                              
			while(fin2>>s)
			{
				if(s=="degrade")
				{
					getline(fin2,s);           
					for(int i = 0;i<ndim;i++)
					{
						getline(fin2,s);
						sscanf (s.c_str(),"%*f %f %f %f %f %f %f %f",
						&temp[0],&temp[1],&temp[2],&temp[3],&temp[4],&temp[5],&temp[6]);
						Wlgth(i) = temp[0]*1000;
						BW(i) = temp[1];
						OPTtrans(i) = temp[2];
						Qeff(i) = temp[3];
						LSB(i) = temp[4];
						FNoise(i) = temp[5];
						Degrade(i) = temp[6];
					}             
				}
			}                                      
		}      
	}

    wNohigh = (int)((10000/(Wlgth(0)/1000-3*BW(0))/5))*5;
    wNolow = (int)((10000/(Wlgth(ndim-1)/1000+3*BW(ndim-1))/5))*5;
    logfile("Sensor construct over");

}

class Processing
{
public:
    Processing(const char* pfilepath);
       
	vector<double> getpfa() const{return pfa;};
	int            getnumfeat() const{return numfeat;};
    int            getatmtype() const{return atmtype;};
    int            getxfrmtype() const{return xfrmtype;};
    int            getalgtype() const{return algtype;};
    int            getmattype() const{return mattype;};
    int            geteiguse() const{return eiguse;};
    vector<double> getos() const{return os;};    
    vector<double> getgrp() const{return grp;};
    int            getbinfact() const{return binfact;};

    void showall(void);	      
       
private:
	vector<double> pfa;
	int            numfeat;
    int            atmtype;
    int            xfrmtype;
    int            algtype;
    int            mattype;
    int            eiguse;    
    vector<double> os;    
    vector<double> grp;
    int            binfact;

};

void Processing::showall(void)
{
     for(size_t i=0;i<pfa.size();i++)
         cout<<"pfa["<<i<<"] = "<<pfa[i]<<endl;
     cout<<"numfeat = "<<numfeat<<endl;
     cout<<"atmtype = "<<atmtype<<endl;
     cout<<"xfrmtype = "<<xfrmtype<<endl;
     cout<<"algtype = "<<algtype<<endl;
     cout<<"mattype = "<<mattype<<endl;
     cout<<"eiguse = "<<eiguse<<endl;
     cout<<"os = "<<os.back()<<endl;
     cout<<"grp = "<<grp.back()<<endl;
     cout<<"binfact = "<<binfact<<endl;                    
}
Processing::Processing(const char* pfilepath)
{
	if(!pfilepath)
	{
		logfile("Processing construct fails, pfilepath NULL","/ERROR");
		throw GeneralException();
	}
	ifstream fin(pfilepath);
	if(!fin.is_open())
	{
		logfile("Processing construct fails, pfilepath Wrong","/ERROR");
		throw GeneralException();
	}
	string s,subs;
	while(fin >> s)
	{
		if(s=="pfa")
		{
			getline(fin,s);
			getline(fin,subs);
			cutcomment(s);
			s = s + subs;
			parsemultiinput(s, pfa);

			fin>>s;
			fin>>s;
			numfeat = atoi(s.c_str()); 
              
 			getline(fin,s);
			fin>>s;
			fin>>s;               
			atmtype = atoi(s.c_str());   
                
 			getline(fin,s);             
			fin>>s;
			fin>>s; 
			xfrmtype = atoi(s.c_str()); 
			    
			getline(fin,s);
 			getline(fin,s);
			fin>>s;
			fin>>s; 
			algtype = atoi(s.c_str());  
                
			getline(fin,s);
			getline(fin,s);
			fin>>s;
			fin>>s;
			mattype = atoi(s.c_str()); 
                
			getline(fin,s);
			getline(fin,s);  
			fin>>s;
			fin>>s;
			eiguse = atoi(s.c_str());
                
			getline(fin,s);  
			getline(fin,s);  

			parsemultiinput(s, os);
              
			getline(fin,s);
			parsemultiinput(s, grp); 
                
			fin>>s;
			fin>>s;
			binfact = atoi(s.c_str());                           
           
			logfile("Processing construct over");

		}

	}

}

void genTape5(Sensor sensor, Scene scene, double surfalb, char* rescode, const char* tape5path, double surftmp = 0)
{
	if(rescode&&tape5path)
	{
		ofstream fout(tape5path);
		if(!fout.is_open())
		{
			logfile("tape5path open wrong","/ERROR");
			throw GeneralException();		
		}
		int imod   = 1, itype  = 2 ,iemsct = 2,
		imult  = 1, im     = 0, noprt  = 0,
		m1     = 0, m2     = 0, m3     = 0,
		m4     = 0, m5     = 0, m6     = 0,
		mdef   = 0, isun   = 0, nstr   = 0,
		iseasn = 0, ivulcn = 0, icstl  = 3,
		irpt   = 0, len    = 0, iph    = 2, 
		ivsa   = 0, isourc = 0, iday   = 180, 
		iparm  = 2, ncralt = -9, ncrspc = -9;
    
    
		double parm1  = 180.0,
		cthik  = -9.0, calt   = -9.0, cext   = -9.0,  //default
		cwavlen= -9.0,  //cloud
		range  = 0,  beta  = 0, ro     = 0,
		ccolwd = -9.0, ccolip = -9.0, chumid = -9.0,  //values
		asymwd = -9.0, asymip = -9.0, co2mx  = 355.0, //hard enter below....     
		rainrt = 0;
    
		char spd    = 'M', dis    = 'F' ,// hard entered below..
		   sun1   = 'F'; 
       
		double tbound = surftmp + 273.15;           
		double angle  = 180.0 - sensor.getviewangle();
		double whh    = scene.getwss();
		double parm2  = scene.getsolangle();
		fout.setf(ios::fixed);
		if (scene.getmodel()==7) im = 1;
		//line 1
		fout<<"TM"<<"  "<<scene.getmodel()<<"    "<<itype<<"    "<<iemsct<<"    "
		  <<imult<<"    "<<m1<<"    "<<m2<<"    "<<m3<<"    "<<m4<<"    " 
		  <<m5<<"    "<<m6<<"    "<<mdef<<"    "<<im<<"    "<<noprt<<
		  setw(8)<<setprecision(3)<<tbound<<
		  setw(8)<<setprecision(2)<<surfalb<<endl;
		//line 2
		fout<<"F   0F   0   355.000                     F T F F"<<endl;
    
		//line 3
		fout<<rescode<<endl;
    
		//line 4
		if(scene.geticld()!=0) scene.seticld(2); 
                  
		fout<<
		setw(5)<<scene.getihaze()<<
		setw(5)<<iseasn<<
		setw(5)<<ivulcn<<
		setw(5)<<icstl<<
		setw(5)<<scene.geticld()<<
		setw(5)<<ivsa<<
		setw(10)<<setprecision(3)<<scene.getmetrange()<<
		setw(10)<<setprecision(3)<<scene.getwss()<<
		setw(10)<<setprecision(3)<<whh<<
		setw(10)<<setprecision(3)<<rainrt<<
		setw(10)<<setprecision(3)<<scene.getgndalt()<<endl;

		if(scene.geticld()!=0)
		{
			fout<<             
			setw(8)<<setprecision(3)<<cthik<<
			setw(8)<<setprecision(3)<<calt<<
			setw(8)<<setprecision(3)<<cext<<
			setw(4)<<ncralt<<
			setw(4)<<ncrspc<<                     
			setw(8)<<setprecision(3)<<cwavlen<<
			setw(8)<<setprecision(3)<<ccolwd<<
			setw(8)<<setprecision(3)<<ccolip<<
			setw(8)<<setprecision(3)<<chumid<<
			setw(8)<<setprecision(3)<<asymwd<<
			setw(8)<<setprecision(3)<<asymip<<endl;                     
		}
	  //line 5
		fout<<
		setw(10)<<setprecision(3)<<sensor.getplatalt()<<
		setw(10)<<setprecision(3)<<scene.getgndalt()<<
		setw(10)<<setprecision(3)<<angle<<
		setw(10)<<setprecision(3)<<range<<
		setw(10)<<setprecision(3)<<beta<<
		setw(10)<<setprecision(3)<<ro<<
		setw(5)<<len<<endl;      
      
	  //line 6
		fout<<
		setw(5)<<iparm<<
		setw(5)<<iph<<
		setw(5)<<iday<<
		setw(5)<<isourc<<endl;
    
	  //line 7
		fout<<
		setw(10)<<setprecision(3)<<parm1<<
		setw(10)<<setprecision(3)<<parm2<<endl; 

	  //line 8     
		fout<<
		setw(10)<<sensor.getwNolow()<<
		setw(10)<<sensor.getwNohigh()<<
		setw(10)<<15<<
		setw(10)<<20<<endl;             
  
	  //line 9
		fout<<
		setw(5)<<irpt<<endl;

		fout.close();
    }
	else
	{
		logfile("genTape5 fails","/ERROR");
		throw GeneralException();
	}
   
}


//ValueB is output at WavelengthB according to WavelengthsA-ValueA, WavelengthA, B both INCREASING
void LinearInterpolation1D(ArrayXd& WavelengthsA, ArrayXd& ValueA, ArrayXd& WavelengthsB, ArrayXd& ValueB)
{
    if(WavelengthsA.size()!=ValueA.size()||WavelengthsB.size()!=ValueB.size()) 
	{
		logfile("LinearInterpolation1D size fails","/ERROR");
		throw GeneralException();
	}
    if(WavelengthsA.tail(1)(0)<WavelengthsB.head(1)(0)||WavelengthsB.tail(1)(0)<WavelengthsB.head(1)(0))
	{
		logfile("LinearInterpolation1D Wlgth increasing fails","/ERROR");
		throw GeneralException();
	}
    if(WavelengthsA.tail(1)(0)<WavelengthsB.tail(1)(0)||WavelengthsA.head(1)(0)>WavelengthsB.head(1)(0))
	{
		logfile("LinearInterpolation1D range fails","/ERROR");
		throw GeneralException();
	}
    int j = 0;   
    for(int i = 0; i < ValueB.size(); i++)
    {
          
        while(WavelengthsB(i)>WavelengthsA(j)) j++;
        ValueB(i) = ValueA(j-1)+(ValueA(j)-ValueA(j-1))/(WavelengthsA(j)-WavelengthsA(j-1))*(WavelengthsB(i)-WavelengthsA(j-1));
    }
    logfile("LinearInterpolation1D Done");

}

//ValueB is output MATRIX at WavelengthB according to WavelengthsA-ValueA
void LinearInterpolation2D(ArrayXd& WavelengthsA, MatrixXd& ValueA, ArrayXd& WavelengthsB, MatrixXd& ValueB)
{
    if(WavelengthsA.size()!=ValueA.cols()||WavelengthsB.size()!=ValueB.rows())
	{
		logfile("LinearInterpolation2D length fails","/ERROR");
		throw GeneralException();
	}
    if(WavelengthsA.tail(1)(0)<WavelengthsB.tail(1)(0)||WavelengthsA.head(1)(0)>WavelengthsB.head(1)(0))
	{
		logfile("LinearInterpolation2D range fails","/ERROR");
		throw GeneralException();
	}

    int m = 0, n = 0;
    double temp1, temp2;   
    for(int i = 0; i < ValueB.rows(); i++)
    {
		while(WavelengthsB(i)>WavelengthsA(m)) m++;
        for(int j = 0; j < i+1; j++)
        {   
            while(WavelengthsB(j)>WavelengthsA(n)) n++;
			////for diagonal element, the interpolation is done by adjancet diagonal elements
            if(m==n&&i==j) ValueB(i,j) = ValueA(m-1,n-1)+(WavelengthsB(j)-WavelengthsA(n-1))*(ValueA(m,n)-ValueA(m-1,n-1))/(WavelengthsA(n)-WavelengthsA(n-1));
            else
            {				
                temp1 = ValueA(m-1,n-1) + (WavelengthsB(j) - WavelengthsA(n-1))*(ValueA(m-1,n) - ValueA(m-1,n-1))/(WavelengthsA(n)-WavelengthsA(n-1));
                temp2 = ValueA(m,n-1) + (WavelengthsB(j) - WavelengthsA(n-1))*(ValueA(m,n) - ValueA(m,n-1))/(WavelengthsA(n)-WavelengthsA(n-1));
                ValueB(i,j) = temp1+(WavelengthsB(i)-WavelengthsA(m-1))*(temp2-temp1)/(WavelengthsA(m)-WavelengthsA(m-1));
                ValueB(j,i) = ValueB(i,j);
            }
//            cout<<m<<" "<<n<<endl;
//            cout<<ValueA(m-1,n-1)<<" "<<ValueA(m-1,n)<<endl;
//            cout<<ValueA(m,n-1)<<" "<<ValueA(m,n)<<endl;            
//            cout<<ValueB(i,j);
//            system("PAUSE");
			
        }
		n = 0;
        
    }

	logfile("LinearInterpolation2D Done");
}

void readWavelengthsMeansandCov(const char* txtpath, ArrayXd& Wavelengths, ArrayXd& Means, MatrixXd& Covs)
{
    if(!txtpath)
	{
		logfile("readWavelengthsMeansandCov txtpath NULL","/ERROR");
		throw GeneralException();
	}
    ifstream fin(txtpath);
    if(!fin.is_open())
	{
		logfile("readWavelengthsMeansandCov txtpath wrong","/ERROR");
		throw GeneralException();
	}
    string s;
    int i = 0;
    for(int j = 0; j<5; j++)
        getline(fin,s);
    fin>>s;fin>>s;fin>>s;fin>>s;
    int sourcedim = atoi(s.c_str());
    Wavelengths.resize(sourcedim);
    Means.resize(sourcedim);
    Covs.resize(sourcedim,sourcedim);
    while( fin >> s ) 
    {  
        if(s == "Wavelengths")
        {
            while(fin >> s && s != "Apparent") 
            {
                Wavelengths(i++) = (double)atof(s.c_str());
            }
   
            while( fin >> s )
            
            {      
                if(s == "Means")
                {                      
                    i = 0;
                    while(fin >> s && s != "Covariances")
                    {
                        double tempMean = (double)atof(s.c_str());      
						if(tempMean>0&&tempMean<=100)
							Means(i++) = (double)atof(s.c_str());
						else if(tempMean<=0)
							Means(i++) = 0.01;
						else if(tempMean>100)
						    Means(i++) = 100;
                    }                                    
                    i = 0;
                    while(fin >> s)
                        Covs(i++) =  (double)atof(s.c_str());

                }                    
                   
            }              
        }     

    }
	logfile("readWavelengthsMeansandCov Done");
}


MatrixXd WriteDiagonalMatrix(double* dataIn, int count)
{
    if(dataIn==NULL||count<1)
	{
		logfile("WriteDiagonalMatrix dataIn NULL","/ERROR");
		throw GeneralException();
	}
    MatrixXd matrixOut(count, count);
    matrixOut.setZero();
    for(int i = 0; i < count ; i++)
    {
        matrixOut(i,i) = dataIn[i];  
    }
    return matrixOut;    
}

MatrixXd WriteDiagonalMatrix(VectorXd dataIn)
{
    int count = dataIn.size();

    MatrixXd matrixOut(count, count);
    matrixOut.setZero();
    for(int i = 0; i < count ; i++)
    {
        matrixOut(i,i) = dataIn(i);  
    }
//    cout<<"WriteDiagonalMatrix"<<endl;
    return matrixOut;    
}

ArrayXd fixedfactor(double f_ifov, double f_aperture,ArrayXd f_transmittance, ArrayXd f_qefficiency, double f_time)
{
      ArrayXd z;
      double PI = 3.1415926f;     
      double HC = 1.986e-16f;
      double A;
      double omega;

      A = PI*(0.5f*f_aperture/10)*(0.5f*f_aperture/10);
	  ////http://en.wikipedia.org/wiki/Solid_angle
      omega = 4*PI*sin(f_ifov/4/1000)*sin(f_ifov/4/1000);
      z = A*omega*f_transmittance*f_qefficiency*f_time/HC;

      return z;
      
}



ArrayXd Sigma_np_s(ArrayXd fixedfactor, ArrayXd L, ArrayXd wlength, ArrayXd bandw)
{
//      fdebug<<"fixedfactor"<<fixedfactor<<endl<<endl<<"L"<<L<<endl<<endl<<"wlength"<<wlength<<endl<<endl<<"bandw"<<bandw<<endl<<endl;  
      ArrayXd output =  fixedfactor*L*wlength*bandw;
      return output;
}

ArrayXd Sigma_nbe_s(int i_nbit, ArrayXd f_LSB, double Be)
{
      int i;
      int length = f_LSB.size();
      ArrayXd z(length);
      z.setZero();
      for (i = 0;i<i_nbit;i++)
	      z = z+(2<<(2*i-1))*f_LSB*f_LSB;
      z = Be/i_nbit*z;
      return z;
}

ArrayXd Sigma_nq_s(ArrayXd f_LSB)
{
      return f_LSB*f_LSB/12;
}

ArrayXd SensorProcessNoise(ArrayXd &SNR, ArrayXd FixedNoise, ArrayXd ffactor, double Cr, double Gn, 
ArrayXd fsigma_nbe_s, ArrayXd fsigma_nq_s,ArrayXd L,ArrayXd WavelengthsSensor,ArrayXd BW, ArrayXd Degrade)
{
	ArrayXd fsigma_np_s = Sigma_np_s(ffactor,L,WavelengthsSensor,BW);
//fdebug<<"fsigma_np_s"<<fsigma_np_s*Degrade<<endl;
	ArrayXd temp = fsigma_np_s*Degrade/L;
	ArrayXd RNoise = (fsigma_np_s*Degrade+FixedNoise*FixedNoise)/(temp*temp);
	ArrayXd CNoise = (L*Cr*L*Cr*0.0001) + fsigma_nbe_s/(temp*temp) + fsigma_nq_s/(temp*temp);
    SNR = fsigma_np_s*Degrade/(sqrt(Gn*Gn*RNoise+CNoise)*temp);
	logfile("SensorProcessNoise done");
	return Gn*Gn*RNoise+CNoise;
}


/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter John Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 */
 
#define LOW 0.02425
#define HIGH 0.97575 
/* Coefficients in rational approximations. */
static const double a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00
};

static const double d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

double invcdf(double p)
{
	double q, r;

	errno = 0;

	if (p < 0 || p > 1)
	{
		errno = EDOM;
		return 0.0;
	}
	else if (p == 0)
	{
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	}
	else if (p == 1)
	{
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	}
	else if (p < LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else if (p > HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}

////bo rewrite
int SearchArray(ArrayXd dataIn, double keyvalue)
{
    for(int i = 0;i<dataIn.size();i++)
    {
        if(dataIn[i]==keyvalue)
           return i;
    }
    return -1;
}

double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}


ArrayXd ROCgentest(ArrayXd Mean, ArrayXd Delta, vector<double> Fraction, vector<double> xposition, int backgroundcount)
{
    int length = 27;
    ArrayXd hm(backgroundcount);
    ArrayXd Pdm(backgroundcount);
    ArrayXd h(27);
    ArrayXi hindex(27);
    ArrayXd Pdmin(27);
    ArrayXd Pfa(27);
    ArrayXd probitP(27);
    int x,y; 
    Pfa.setZero();
    ArrayXd xpositiongrid(27);
    xpositiongrid<<1e-8, 2e-8, 5e-8,
                   1e-7, 2e-7, 5e-7,
                   1e-6, 2e-6, 5e-6,
                   1e-5, 2e-5, 5e-5,
                   1e-4, 2e-4, 5e-4,
                   1e-3, 2e-3, 5e-3,
                   1e-2, 2e-2, 5e-2,
                   1e-1, 2e-1, 5e-1,
                   1,2,5;
//1e-9, 2e-9, 5e-9,
    for(int i = 0; i< 27; i++)
    {
        probitP[i] = invcdf(1-xpositiongrid(length-1-i));
//cout<<"ProbitP="<<probitP[i]<<endl ;
    }
         
    Mean[backgroundcount] = 0.292808;
    Delta[backgroundcount] = 0.0554;
    Mean[0] = -0.0027;
    Mean[1] = -0.015689;
    Mean[2] = 0.07628;
    Delta[0] = 0.0135 ;   
    Delta[1] = 0.0492 ;   
    Delta[2] = 0.0428;
    				
    Fraction[0] *=0.01;
    Fraction[1] *=0.01;
    Fraction[2] *=0.01;	
//cout<<"Fraction[0]"<<Fraction[0]<<endl;	
//cout<<"Fraction[1]"<<Fraction[1]<<endl;	
//cout<<"Fraction[2]"<<Fraction[2]<<endl;		
    for (y=0;y<27;y++)
    {
        for (x=0;x<backgroundcount;x++)
        {
            hm[x] = Mean[x]+Delta[x]*probitP[y]	;	
            Pdm[x] = 0.5 - 0.5*erf((hm[x]-Mean[backgroundcount])/(1.414*Delta[backgroundcount]));
        }

        hindex[y] = SearchArray(Pdm, Pdm.minCoeff());

//      hindex[y] = Pdm.index(Pdm.minCoeff())
        h[y] = hm[hindex[y]];
        Pdmin[y] = Pdm[hindex[y]];        
        
        for(x=0;x<backgroundcount;x++)
            Pfa[y] = Pfa[y] + Fraction[x]*(0.5 - 0.5*erf((h[y]-Mean[x])/(1.414*Delta[x]))) ;
//cout<<"["<<x+y*3<<"]="<<0.5 - 0.5*erf((h[y]-Mean[x])/(1.414*Delta[x]))<<endl;            }
    }
    
//    cout<<erf(0)<<" "<<erf(0.5)<<""<<erf(1)<<""<<erf(-0.5)<<""<<erf(-1)<<""<<endl; 
//    fROC<<Pdmin<<endl;
//    fROC<<Pfa<<endl;
    return Pfa;     
        
}

ArrayXd ROCgen(ArrayXd Mean, ArrayXd Delta, vector<double> Fraction, vector<double> xposition, int backgroundcount)
{
    int length = 72; //for pfa grid calculation
	int xlength = xposition.size();
    ArrayXd hm(backgroundcount);
    ArrayXd Pdm(backgroundcount);

    ArrayXd xpositiongrid(length);
    xpositiongrid<<1.5e-8, 2.5e-8, 3.5e-8, 4.5e-8, 5.5e-8, 6.5e-8, 7.5e-8, 8.5e-8, 9.5e-8,
                   1.5e-7, 2.5e-7, 3.5e-7, 4.5e-7, 5.5e-7, 6.5e-7, 7.5e-7, 8.5e-7, 9.5e-7,
                   1.5e-6, 2.5e-6, 3.5e-6, 4.5e-6, 5.5e-6, 6.5e-8, 7.5e-8, 8.5e-6, 9.5e-6,
                   1.5e-5, 2.5e-5, 3.5e-5, 4.5e-5, 5.5e-5, 6.5e-5, 7.5e-5, 8.5e-5, 9.5e-5,
                   1.5e-4, 2.5e-4, 3.5e-4, 4.5e-4, 5.5e-4, 6.5e-4, 7.5e-4, 8.5e-4, 9.5e-4,
                   1.5e-3, 2.5e-3, 3.5e-3, 4.5e-3, 5.5e-3, 6.5e-3, 7.5e-3, 8.5e-3, 9.5e-3,
                   1.5e-2, 2.5e-2, 3.5e-2, 4.5e-2, 5.5e-2, 6.5e-2, 7.5e-2, 8.5e-2, 9.5e-2,
                   1.5e-1, 2.5e-1, 3.5e-1, 4.5e-1, 5.5e-1, 6.5e-1, 7.5e-1, 8.5e-1, 9.5e-1;

    ArrayXd h(length);
    ArrayXi hindex(length);
    ArrayXd Pdmin(length);
    ArrayXd Pfa(length);
    ArrayXd probitP(length);
	ArrayXd xpositionarray(xlength);
	ArrayXd outPfa(xlength);
	ArrayXd outPdmin(xlength);
	ArrayXd outPdPfa = ArrayXd::Zero(length*2+2);
    int x,y; 
    Pfa.setZero();

    for(x=0;x<backgroundcount;x++)
        Fraction[x] *=0.01;
	for(x=0;x<xlength;x++)
		xpositionarray(x) = xposition[xlength-x-1];

    for(int i = 0; i< length; i++)
        probitP[i] = invcdf(1-xpositiongrid[length-1-i]);
/*         
    Mean[backgroundcount] = 0.292808;
    Delta[backgroundcount] = 0.0554;
    Mean[0] = -0.0027;
    Mean[1] = -0.015689;
    Mean[2] = 0.07628;
    Delta[0] = 0.0135 ;   
    Delta[1] = 0.0492 ;   
    Delta[2] = 0.0428;
*/  	
    for (y=0;y<length;y++)
    {
        for (x=0;x<backgroundcount;x++)
        {
			if(Mean[x] > Mean[backgroundcount]) 
				Mean[x] = 2 * Mean[backgroundcount] - Mean[x];
            hm[x] = Mean[x]+Delta[x]*probitP[y];
			
            Pdm[x] = 0.5 - 0.5*erf((hm[x]-Mean[backgroundcount])/(1.414*Delta[backgroundcount]));
        }
//		cout<<"hm["<<y<<"] = "<<endl<<hm<<endl;
//		cout<<"Pdm["<<y<<"] = "<<endl<<Pdm<<endl;
        hindex[y] = SearchArray(Pdm, Pdm.minCoeff());
//         hindex[y] = Pdm.index(Pdm.minCoeff())
        h[y] = hm[hindex[y]];
        Pdmin[y] = Pdm[hindex[y]];
        for(x=0;x<backgroundcount;x++)
        {
            Pfa[y] = Pfa[y] + Fraction[x]*(0.5 - 0.5*erf((h[y]-Mean[x])/(1.414*Delta[x]))) ;
        }
    } 
    LinearInterpolation1D(xpositiongrid, Pfa, xpositionarray, outPfa);
	LinearInterpolation1D(xpositiongrid, Pdmin, xpositionarray, outPdmin);
	outPdPfa<<1,Pdmin,1,Pfa;
    return outPdPfa;
        
}

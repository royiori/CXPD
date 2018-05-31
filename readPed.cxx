//____________________________________________
// 1个数据文件包含8个channel
// 1个channel包含72x72个pixel的数据
// 1个pixel的数据包含从第0frame到第808frame的数据顺序排列
// so:
// 1个数据文件大小为: 72x72x809
//
// 采集回来的.pd1文件中只有一个channel是有数据的
// 将其他不用的数据去掉，并按frame重新排列
// 
//  .pd1   ---->   .mdat
//
//                               - lq
//____________________________________________


const int nFrame = 808;
const int nChannel = 8;
const int nXpixel = 72;
const int nYpixel = 72;
const int nLength = (nXpixel*nYpixel)*nFrame;
const int channel = 0;
unsigned short data[nLength];

const char *filetypes[] = {
                            "raw data file", "*.pd1",
                            0,               0 };


//-----
// 
TString SelectFile()
{
  static TString dir(".");
  TGFileInfo fi;
  fi.fFileTypes = filetypes;
  fi.fIniDir    = StrDup(dir);
  printf("fIniDir = %s\n", fi.fIniDir);
  new TGFileDialog(gClient->GetRoot(), NULL, kFDOpen, &fi);
  if(fi.fFilename==NULL) { printf("No file selected!\n"); return TString("NONE"); }
  printf("Open directory: %s\n", fi.fIniDir);
  return fi.fIniDir;
}

//___________________________
// main
void readPed()
{
    TString path = SelectFile();
    FILE *fp = gSystem->OpenPipe("ls "+path+"/*.pd1", "r");
    if(!fp) { cout<<"----> NO pd1 data exists in "<<path<<"!"<<endl; return;}
    
    vector <TString> pdList;    
    char line[1000];
    while(fgets(line, sizeof(line), fp))
    {
        TString s(line);
        if(s.Index(".pd1")==-1) continue;
        pdList.push_back(s.ReplaceAll("\n",""));
    }
    cout<<"----> "<<pdList.size()<<" pd1 data files exist:"<<endl;


    for(int i=0; i<(int)pdList.size(); i++)
    {
        TString mdatName = pdList[i];
        mdatName.Replace(mdatName.Index(".pd1"), 4, ".mdat");

        if(!gSystem->AccessPathName(mdatName, kFileExists)) 
        {
            cout<<"-----> "<<mdatName<<" is existed."<<endl;
            continue;
        }
        cout<<"----> converting "<<pdList[i]<<" to "<<mdatName<<endl;;    
        
        ifstream ifSignal(pdList[i], ios::binary);

        ofstream oftest(mdatName, ios::binary);

        unsigned short _data;
        //vector<vector <int>> bgdata;
        //bgdata.resize(nFrame);
        //for(int i=0; i<nFrame; i++) bgdata[i].resize(nXpixel*nYpixel);

        if(ifSignal.good())
        {
            // read data
            //
            ifSignal.seekg( (channel - nChannel) * nLength * sizeof(unsigned short), ios::end);
            ifSignal.read((char *)data, sizeof(data));

            for(int ii = 0; ii < nFrame; ii++)
            {
                for(int jj = 0; jj<nXpixel*nYpixel; jj++)
                {
                    oftest.write((char *)(&data[ jj * nFrame + ii ]), sizeof(unsigned short));
                }
            }

            /*for(int ii=0; ii<nFrame; ii++) 
	        {
                for(int jj=0; jj<(nXpixel*nYpixel); jj++) 
	            {
                    ifSignal.seekg((-(nChannel*(nXpixel*nYpixel*(nFrame-ii)-jj))+channel)*sizeof(_data),ios::end);
                    ifSignal.read((char*)(&_data), sizeof(_data));
                    //bgdata[i][j] = _data;
		            oftest.write((char *)(&_data), sizeof(_data));
                }
            }*/
        }

        ifSignal.close();
        oftest.close();
    }
}






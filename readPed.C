// dat[0] = 65535                  ;;  required flag
// dat[1] = xmin & dat[2] = xmax    ;;  x min and max in pixel
// dat[3] = ymin & dat[4] = ymax     ;;  y min and max in pixel
// dat[5] = evtid                            ;; event id
// dat[6] = 0 & dat[7] = 0               ;;  time stamp
// dat[8] = 0 & dat[9] = 0               ;;   time stamp
// dat[10:*] = ADC                       ;; ADC

void readPed()
{
    ifstream ifSignal("./data/out835.pd1", ios::binary);

    ofstream oftest("./data/out835.mdat", ios::binary);


    const int nFrame = 808;
    const int nChannel = 8;
    const int nXpixel = 72;
    const int nYpixel = 72;
    const int nLength = nChannel*(nXpixel*nYpixel)*nFrame;
    const int nHeader = 1024;

    int headerdat[nHeader];
    vector<vector <int>> bgdata;
    unsigned short _data;
    int channel = 0;
 
    bgdata.resize(nFrame);
    for(int i=0; i<nFrame; i++) bgdata[i].resize(nXpixel*nYpixel);

    if(ifSignal.good())
    {
        for(int i=0; i<nHeader; i++) 
	{
            ifSignal.read((char*)(&_data), sizeof(_data));
            headerdat[i] = _data;
        }

        //
        // readout the background data
        //
        for(int i=0; i<nFrame; i++) 
	{
            for(int j=0; j<(nXpixel*nYpixel); j++) 
	    {
                ifSignal.seekg((-(nChannel*(nXpixel*nYpixel*(nFrame-i)-j))+channel)*sizeof(_data),ios::end);
                ifSignal.read((char*)(&_data), sizeof(_data));
                bgdata[i][j] = _data;
		oftest.write((char *)(&_data), sizeof(_data));
            }
        }
    }

    ifSignal.close();
    oftest.close();
}





////------------
////PIXE:
//    unsigned short _data;
//    short xmin, xmax, ymin, ymax;
//
//
//
//    for (int ii = 0; ii < 10; ii++)
//    {
//        for (int j = 0; j <= 9; j++)
//        {
//            ifSignal.read((char *)(&_data), sizeof(_data));
//            if (j == 1) 
//            {
//                ymin = _data;
//                _data = 0;
//                //cout<<"ymin: "<<ymin<<"->"<<_data<<endl;
//            }
//            if (j == 2)
//            {
//                ymax = _data;
//                _data -= (ymin-0);
//                //cout<<"ymax: "<<ymax<<"->"<<_data<<endl;
//            }
//            if (j == 3)
//            {
//                xmin = _data;
//                _data = 0;
//                //cout<<"xmin: "<<xmin<<"->"<<_data<<endl;
//            }
//            if (j == 4)
//            {
//                xmax = _data;
//                _data -= (xmin-0);
//                //cout<<"xmax: "<<xmax<<"->"<<_data<<endl;
//            }
//            oftest.write((char *)(&_data), sizeof(_data));
//        }
//
//        for (int i = 0; i <= xmax-xmin; i++)
//        {
//            for (int j = 0; j <= ymax-ymin; j++)
//            {
//                ifSignal.read((char *)(&_data), sizeof(_data));
//                if(_data<240) _data = 0;
//                else {cout<<ii<<"-- "<<i<<", "<<j<<" : "<<_data<<endl ; _data = 100;}
//                //else if(_data>=8 && _data<20) _data = 30;
//                //else if(_data>=20)            _data = 50;
//                oftest.write((char *)(&_data), sizeof(_data));
//            }
//        }
//    }
//    /*
//    for (int ii = 0; ii < 2; ii++)
//    {
//        data = 65535;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 30;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 40;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 30;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 40;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 0;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 0x2fee;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 0;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 0x2fee;
//        oftest.write((char *)(&data), sizeof(data));
//        data = 0x2fee;
//        oftest.write((char *)(&data), sizeof(data));
//        for (int i = 30; i < 40; i++)
//        {
//            for (int j = 30; j < 40; j++)
//            {
//                data = 100;
//                oftest.write((char *)(&data), sizeof(data));
//            }
//        }
//    }
//    return;
//    */


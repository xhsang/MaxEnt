PVOID load_image(char *Filename,int &x,int &y)
{
	HANDLE FileToOpen=::CreateFileA(Filename,GENERIC_READ,FILE_SHARE_READ,
		NULL,OPEN_ALWAYS,FILE_ATTRIBUTE_NORMAL,NULL); 
	if(FileToOpen==NULL)
	{
		printf("loading bitmap file error!\n");
		return NULL;
	}
	DWORD filesize,filesizehigh;
	filesize=GetFileSize(FileToOpen,&filesizehigh);
	
	PVOID p;
	p=VirtualAlloc(NULL,filesize,MEM_COMMIT,PAGE_READWRITE);
	if(p==NULL)
	{
		printf("alloc memory for image error!\n");
		CloseHandle(FileToOpen);
		return NULL;
	}
	if(!ReadFile(FileToOpen,p,filesize,&filesizehigh,NULL))
	{
		printf("Read imamge file error!\n");
		CloseHandle(FileToOpen);
		VirtualFree(p,filesize,MEM_DECOMMIT);
		return NULL;
	}
	CloseHandle(FileToOpen);

	BITMAPFILEHEADER *pFileHeader;
	pFileHeader=(BITMAPFILEHEADER*)p;
	BITMAPINFOHEADER *pInfoHeader;
	pInfoHeader=(BITMAPINFOHEADER*)((BYTE*)p+sizeof(BITMAPFILEHEADER));
	//DWORD d1=sizeof(BITMAPINFOHEADER );  //0x28
	//DWORD d2=sizeof(BITMAPV4HEADER);     //0x6c
	//DWORD d3=sizeof(BITMAPV5HEADER);     //0x7c
	if(pInfoHeader->biSize!=sizeof(BITMAPINFOHEADER))
	{
		printf("Cannot load image!Please make sure that the image is of 24-bits bitmap");
		VirtualFree(p,filesize,MEM_DECOMMIT);
		return NULL;
	}
	x=pInfoHeader->biWidth;
	y=pInfoHeader->biHeight;
	return p;
}

int set_image_data_to_b(PVOID p,float *b,float *sigma)
{
	BITMAPFILEHEADER *pFileHeader;
	pFileHeader=(BITMAPFILEHEADER*)p;
	BITMAPINFOHEADER *pInfoHeader;
	pInfoHeader=(BITMAPINFOHEADER*)((BYTE*)p+sizeof(BITMAPFILEHEADER));
	BYTE *pdata=(BYTE*)((BYTE*)p+pFileHeader->bfOffBits);
	int x,y,i,j;
	x=pInfoHeader->biWidth;
	y=pInfoHeader->biHeight;
	for(i=0;i<x;i++)
	{
		for(j=0;j<y;j++)
		{
			*(b+i*y+j)=(float)(*(pdata+3*(i*y+j)));
			BYTE t1=(*(pdata+3*(i*y+j)));
			BYTE t2=(*(pdata+3*(i*y+j)+1));
			BYTE t3=(*(pdata+3*(i*y+j)+2));
			if(t1!=t2||t1!=t3||t2!=t3)
			{
				printf("not a gray bimap,cannot load\n");
				return 0;
			}
			*(sigma+i*y+j)=5;
		}
	}
	return pFileHeader->bfOffBits;
}

int save_image_data(PVOID pimage,int times,float alfa,char *directory)
{
	char filename[100];
	sprintf(filename,"%s\\time%d_alfa%f.bmp",directory,times,alfa);
	HANDLE FileToOpen=CreateFileA(filename,GENERIC_WRITE,0,NULL,CREATE_ALWAYS,0,NULL);
	DWORD reserve;
	BITMAPFILEHEADER *pFileHeader;
	pFileHeader=(BITMAPFILEHEADER*)pimage;
	WriteFile(FileToOpen,pimage,pFileHeader->bfSize,&reserve,NULL);
	CloseHandle(FileToOpen);
	return 0;
}

int set_gray_pallette(PVOID p)
{
	BITMAPFILEHEADER *pFileHeader;
	pFileHeader=(BITMAPFILEHEADER*)p;
	BITMAPINFOHEADER *pInfoHeader;
	pInfoHeader=(BITMAPINFOHEADER*)((BYTE*)p+sizeof(BITMAPFILEHEADER));
	BYTE *pdata=(BYTE*)((BYTE*)p+pFileHeader->bfOffBits);
	BYTE *pallette=(BYTE*)((BYTE*)p+sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER));
	int i,j,color;
	BYTE R,G,B,Z;
	for(i=0;i<pInfoHeader->biWidth;i++)
	{
		for(j=0;j<pInfoHeader->biHeight;j++)
		{
			color=*(pdata+i*pInfoHeader->biHeight+j);
			R=*(pallette+color*4);
			G=*(pallette+color*4+1);
			B=*(pallette+color*4+2);
			Z=*(pallette+color*4+3);
			*(pdata+i*pInfoHeader->biHeight+j)=BYTE((R+G+B)/3.0);
		}
	}
	for(i=0;i<256;i++)
	{
		*(pallette+i*4)=i;
		*(pallette+i*4+1)=i;
		*(pallette+i*4+2)=i;
		*(pallette+i*4+3)=0;
	}
	return 0;
}

int create_directory_for_images(char *directory,char *p_directory,int mx,int my,float sigma_x,float sigma_y)
{
	DWORD t=GetTickCount();
	//sprintf(directory,"%s\\%d_probe%d-%dHW%2.1fp_t_d%2.1f",p_directory,t,mx,my,sigma_x,sigma_y);
	sprintf(directory,"%s\\%d",p_directory,t);
	if(CreateDirectoryA(directory,NULL)==0)
	{
		DWORD err=GetLastError();
		printf("create directory error!\n");
		return 0;
	}
	return 1;
}
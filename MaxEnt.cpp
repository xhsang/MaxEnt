// MaxEnt.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "math.h"
#include "stdio.h"
#include "windows.h"
#include "stdlib.h"
#include "handle_bitmap.h"

#define b(i,j) (*(b+((int)(i))*N[1]+(int)(j)))
#define sigma(i,j) (*(sigma+((int)(i))*N[1]+(int)(j)))
#define h(i,j) (*(h+((int)(i))*n[1]+(int)(j)))
#define h0(i,j) (*(h0+((int)(i))*n[1]+(int)(j)))
#define dh(i,j) (*(dh+((int)(i))*n[1]+(int)(j)))
#define g(i,j) (*(g+((int)(i))*n[1]+(int)(j)))
#define mu(i,j) (*(mu+((int)(i))*n[1]+(int)(j)))
#define c(i,j) (*(c+((int)(i))*m[1]+(int)(j)))
#define A(i,j,k) (*(A+(i)*n[1]*(2*m[0]-1)*(2*m[1]-1)+(j)*(2*m[0]-1)*(2*m[1]-1)+k))
#define B(i,j,k) (*(B+(i)*n[1]*(2*m[0]-1)*(2*m[1]-1)+(j)*(2*m[0]-1)*(2*m[1]-1)+k))
#define Max_Iteration 500
#define ERR 0.001

int n[2],m[2],N[2];   //n:sizeof(h); m:sizeof(probe); N:sizeof(picture). Obviously, n=m+N-1;
float omiga=0.12;      //SOR factor
struct PSF
{
	float Cs;
	float lamda;
	float alfa;
	float df;
	float gamma; //gamma equates to width of half intensity
	float pixel_to_distance;
}psf;

int alloc_initiate(float *&h,float *&h0,float *&dh,float *&mu,float *&g,float *&c,float *&b,float *&A,float *&B,float *&sigma)
{
	h=(float*)VirtualAlloc(NULL,sizeof(float)*n[0]*n[1],MEM_COMMIT,PAGE_READWRITE);
	if(h==NULL)
		return 0;
	h0=(float*)VirtualAlloc(NULL,sizeof(float)*n[0]*n[1],MEM_COMMIT,PAGE_READWRITE);
	if(h0==NULL)
		return 0;
	dh=(float*)VirtualAlloc(NULL,sizeof(float)*n[0]*n[1],MEM_COMMIT,PAGE_READWRITE);
	if(dh==NULL)
		return 0;
	mu=(float*)VirtualAlloc(NULL,sizeof(float)*n[0]*n[1],MEM_COMMIT,PAGE_READWRITE);
	if(mu==NULL)
		return 0;
	g=(float*)VirtualAlloc(NULL,sizeof(float)*n[0]*n[1],MEM_COMMIT,PAGE_READWRITE);
	if(g==NULL)
		return 0;
	c=(float*)VirtualAlloc(NULL,sizeof(float)*m[0]*m[1],MEM_COMMIT,PAGE_READWRITE);
	if(c==NULL)
		return 0;
	b=(float*)VirtualAlloc(NULL,sizeof(float)*N[0]*N[1],MEM_COMMIT,PAGE_READWRITE);
	if(b==NULL)
		return 0;
	sigma=(float*)VirtualAlloc(NULL,sizeof(float)*N[0]*N[1],MEM_COMMIT,PAGE_READWRITE);
	if(sigma==NULL)
		return 0;
	A=(float*)VirtualAlloc(NULL,sizeof(float)*n[0]*n[1]*(2*m[0]-1)*(2*m[1]-1),MEM_COMMIT,PAGE_READWRITE);
	if(A==NULL)
		return 0;
	//B=(float*)VirtualAlloc(NULL,sizeof(float)*n[0]*n[1]*(2*m[0]-1)*(2*m[1]-1),MEM_COMMIT,PAGE_READWRITE);
	B=(float*)VirtualAlloc(NULL,sizeof(float)*n[0]*n[1],MEM_COMMIT,PAGE_READWRITE);
	if(B==NULL)
		return 0;
	return 1;
}

int alloc_free(float *h,float *h0,float *dh,float *mu,float*g,float *c,float *b,float *A,float *B,float *sigma)
{
	if(VirtualFree((PVOID)h,sizeof(float)*n[0]*n[1],MEM_DECOMMIT)==NULL)
		return 0;
	if(VirtualFree((PVOID)h0,sizeof(float)*n[0]*n[1],MEM_DECOMMIT)==NULL)
		return 0;
	if(VirtualFree((PVOID)dh,sizeof(float)*n[0]*n[1],MEM_DECOMMIT)==NULL)
		return 0;
	if(VirtualFree((PVOID)mu,sizeof(float)*n[0]*n[1],MEM_DECOMMIT)==NULL)
		return 0;
	if(VirtualFree((PVOID)g,sizeof(float)*n[0]*n[1],MEM_DECOMMIT)==NULL)
		return 0;
	if(VirtualFree((PVOID)c,sizeof(float)*m[0]*m[1],MEM_DECOMMIT)==NULL)
		return 0;
	if(VirtualFree((PVOID)b,sizeof(float)*N[0]*N[1],MEM_DECOMMIT)==NULL)
		return 0;
	if(VirtualFree((PVOID)sigma,sizeof(float)*N[0]*N[1],MEM_DECOMMIT)==NULL)
		return 0;
	if(VirtualFree((PVOID)A,sizeof(float)*n[0]*n[1]*(2*m[0]-1)*(2*m[1]-1),MEM_DECOMMIT)==NULL)
		return 0;
	//if(VirtualFree((PVOID)B,sizeof(float)*n[0]*n[1]*(2*m[0]-1)*(2*m[1]-1),MEM_DECOMMIT)==NULL)
	//	return 0;
	if(VirtualFree((PVOID)B,sizeof(float)*n[0]*n[1],MEM_DECOMMIT)==NULL)
		return 0;
	return 1;
}


int print_matrix(float *p,int l1,int l2)
{
	for(int i=0;i<l1;i++)
	{
		for(int j=0;j<l2;j++)
			printf("%7.6f  ",*(p+i*l1+j));
		printf("\n");
	}
	printf("\n\n");
	return 0;
}

float calculate_PSF(float Cs,float lamda,float alfa,float df,float r0,float sita0)
{
	//calculate probe spread function in a specific positon(r0,sita0)
	//use simpson formula
	int m,n,co,i,j;
	float dR,dsita,Rmax,integral_real,integral_image,gamma,inner,R,sita;
	float pi=3.1415926;
	m=50;
	n=50;
	Rmax=sin(alfa)/lamda;
	dR=Rmax/(2*m);
	dsita=(2*pi)/(2*n);
	integral_real=integral_image=0;
	for(R=0,i=0;R<=Rmax;R+=dR,i++)
	{
		for(sita=0,j=0;sita<=2*pi;sita+=dsita,j++)
		{
			if((i==0||i==2*m)&&(j==0||j==2*n))
				co=1;
			else if((i==0)||(i==2*m))
			{
				if(j%2==1)
					co=4;
				else
					co=2;
			}
			else if((j==0)||(j==2*n))
			{
				if(i%2==1)
					co=4;
				else
					co=2;
			}
			else
			{
				co=4;
				co*=(i%2)+1;
				co*=(j%2)+1;
			}
			gamma=pi*lamda*R*R*(df+(Cs*lamda*lamda*R*R)/2);
			inner=2*pi*r0*R*cos(sita-sita0);
			integral_real+=co*R*cos(inner+gamma);
			integral_image+=co*R*sin(inner+gamma);
		}
	}
	integral_real*=dR*dsita/9;
	integral_image*=dR*dsita/9;
	return (integral_real*integral_real+integral_image*integral_image);
}

int set_probe(float *c,float radius)
{
	int i,j;
	for(i=0;i<m[0];i++)
		for(j=0;j<m[1];j++)
			c(i,j)=((m[0]-1)/2.+(m[1]-1)/2.-abs(i-(m[0]-1)/2.)-abs(j-(m[1]-1)/2.))*radius;
	//print_matrix(c,m[0],m[1]);
	return 1;
}

int set_normal_dist_probe(float *c,float sigma_x,float sigma_y)
{
	if(sigma_x==0&&sigma_y==0)
	{
		sigma_x=9;
		sigma_y=9;
	}
	int i,j,s,t;
	for(i=0;i<m[0];i++)
	{
		for(j=0;j<m[1];j++)
		{
			s=i-(m[0]-1)/2;
			t=j-(m[1]-1)/2;
			c(i,j)=exp(-(s*s/(sigma_x*sigma_x)+t*t/(sigma_y*sigma_y)));
		}
	}
	//print_matrix(c,m[0],m[1]);
	return 1;
}

int set_PSF_probe(float *c)
{
	int i,j,s,t;
	float r0,sita0;
	for(i=0;i<m[0];i++)
	{
		for(j=0;j<m[1];j++)
		{
			s=i-(m[0]-1)/2;
			t=j-(m[1]-1)/2;
			sita0=0;
			r0=sqrt((float)(s*s+t*t))/psf.pixel_to_distance;
			c(i,j)=calculate_PSF(psf.Cs,psf.lamda,psf.alfa,psf.df,r0,sita0);
		}
	}
	print_matrix(c,m[0],m[1]);
	return 1;
}

int set_Lorenze_probe(float *c)
{
	int i,j,s,t;
	float r0;
	for(i=0;i<m[0];i++)
	{
		for(j=0;j<m[1];j++)
		{
			s=i-(m[0]-1)/2;
			t=j-(m[1]-1)/2;
			r0=sqrt((float)(s*s+t*t))/psf.pixel_to_distance;
			c(i,j)=psf.gamma/(3.1415926*(r0*r0+psf.gamma*psf.gamma));
		}
	}
	print_matrix(c,m[0],m[1]);
	return 0;
}
int get_image(float *b,float *h,float *c,float *sigma)   
{
	int i,j;
	for(i=0;i<n[0];i++)
		for(j=0;j<n[1];j++)
			h(i,j)=0;

	//initialize picture to be deconvoluted,using self_constructed picture first
	//for(i=n[0]/2+m[0];i<n[0]-m[0];i+=2*m[0])
	//{
	//	for(j=n[1]/2+m[1];j<n[1]-m[1];j+=2*m[1])
	//	{
	//		h(i,j)=255;
	//		h(i,n[1]-j)=255;
	//		h(n[0]-i,j)=255;
	//		h(n[0]-i,n[1]-j)=255;
	//	}
	//}
	int number=20,i1,j1,range=4,overlap;
	srand(GetTickCount());
	do{
		overlap=0;
		i=rand()%n[0];
		j=rand()%n[1];
		if(i>(n[0]-(m[0]-1))||i<(m[0]-1))
			continue;
		if(j>(n[1]-(m[1]-1))||j<(m[1]-1))
			continue;
		for(i1=-range;i1<=range;i1++)
		{
			for(j1=-range;j1<=range;j1++)
			{
				if(h(i+i1,j+j1)==255)
				{
					overlap=1;
					break;
				}
			}
			if(overlap==1)
				break;
		}
		if(overlap==1)
			continue;
		h(i,j)=255;
		number--;
	}while(number>0);
	//print_matrix(h,n[0],n[1]);

	set_normal_dist_probe(c,0,0);

	
	for(i=0;i<N[0];i++)
	{
		for (j=0;j<N[1];j++)
		{
			int i1,j1;
			b(i,j)=0;
			for(i1=0;i1<m[0];i1++)
				for(j1=0;j1<m[1];j1++)
					b(i,j)+=h(i+i1,j+j1)*c(m[0]-1-i1,m[1]-1-j1);
			sigma(i,j)=1;
		}
	}
	//print_matrix(b,N[0],N[1]);
	return 1;
}

int set_starting_point(float *h0)
{
	int i,j;
	for(i=0;i<n[0];i++)
		for(j=0;j<n[1];j++)
			h0(i,j)=1.0;
	//print_matrix(h0,n[0],n[1]);
	return 1;
}

int get_g(float *g,float *h,float *h0,float *c,float *b,float *sigma,float alfa)
{
	int s,t,p,q,u,v;
	float x;
	for(s=0;s<n[0];s++)
	{
		for(t=0;t<n[1];t++)
		{
			g(s,t)=0;
			for(p=s+1-m[0];p<=s;p++)
			{
				if(p<0||p>N[0]-1)
					continue;
				for(q=t+1-m[1];q<=t;q++)
				{
					if(q<0||q>N[1]-1)
						continue;
					x=0;
					for(u=p;u<=m[0]+p-1;u++)
					{
						if(u<0||u>n[0]-1)
							continue;
						for(v=q;v<=m[1]+q-1;v++)
						{
							if(v<0||v>n[1]-1)
								continue;
							x+=h(u,v)*c(m[0]-u+p-1,m[1]-v+q-1);
						}
					}
					g(s,t)-=c(m[0]-s+p-1,m[1]-t+q-1)*(x-b(p,q))/(sigma(p,q)*sigma(p,q));
				}
			}
			g(s,t)-=alfa*log(h(s,t)/h0(s,t));
		}
	}
	//print_matrix(g,n[0],n[1]);
	return 1;
}

int get_L_gradient2(float *A,float *c)
{
	int s,t,x,y,p,q,count;
	for(s=0;s<n[0];s++)
	{
		for(t=0;t<n[1];t++)
		{
			count=0;
			for(x=s+1-m[0];x<=s+m[0]-1;x++)
			{
				if(x<0||x>n[0]-1)
					continue;
				for(y=t+1-m[1];y<=t+m[1]-1;y++)
				{
					if(y<0||y>n[1]-1)
						continue;
					A(s,t,count)=0;
					for(p=max(s+1-m[0],x+1-m[0]);p<=min(s,x);p++)
					{
						if(p<0||p>N[0]-1)
							continue;
						for(q=max(t+1-m[1],y+1-m[1]);q<=min(t,y);q++)
						{
							if(q<0||q>N[1]-1)
								continue;
							A(s,t,count)+=c(m[0]-s+p-1,m[1]-t+q-1)*c(m[0]-x+p-1,m[1]-y+q-1);
						}
					}
					count++;
				}
			}
		}
	}
	return 1;
}

int get_AB(float *A,float *B,float *c,float *mu,float beta)
{
	int s,t,x,y,count;
	for(s=0;s<n[0];s++)
	{
		for(t=0;t<n[1];t++)
		{
			count=0;
			for(x=s+1-m[0];x<=s+m[0]-1;x++)
			{
				if(x<0||x>n[0]-1)
					continue;
				for(y=t+1-m[1];y<=t+m[1]-1;y++)
				{
					if(y<0||y>n[1]-1)
						continue;
					B(s,t,count)=A(s,t,count)*sqrt(mu(s,t)*mu(x,y));
					if(s==x&&t==y)
						B(s,t,count)+=beta;
					count++;
				}
			}
		}
	}
	return 1;
}

int get_dh_SOR(float *B,float *mu,float *g,float *dh,float *h)
{
	int i,j;
	for(i=0;i<n[0];i++)
	{
		for(j=0;j<n[1];j++)
		{
			if(mu(i,j)<0)
			{
				printf("error!sqrt of a negative number!\n");
				return 0;
			}
			g(i,j)*=sqrt(mu(i,j));
			dh(i,j)=1;           //designating dh initial value
		}
	}
	int s,t,count,iteration;
	float divider,temp,err;
	iteration=1;
	do{
		err=0;
		for(i=0;i<n[0];i++)
		{
			for(j=0;j<n[1];j++)
			{
				count=0;
				divider=0;
				temp=dh(i,j);
				dh(i,j)=g(i,j);
				for(s=max(0,i+1-m[0]);s<=min(n[0]-1,i+m[0]-1);s++)
				{
					for(t=max(0,j+1-m[1]);t<=min(n[1]-1,j+m[1]-1);t++)
					{
						if(i==s&&j==t)
						{
							divider=B(i,j,count);
							dh(i,j)-=B(i,j,count)*temp;
						}
						else
							dh(i,j)-=B(i,j,count)*dh(s,t);
						count++;
					}
				}
				if(divider==0)
					printf("error when calculating B!\n");
				dh(i,j)/=divider;
				dh(i,j)*=omiga;
				dh(i,j)+=temp;
				temp=abs(temp-dh(i,j));
				if(temp>err)
					err=temp;
			}
		}
		iteration++;
	}while(iteration<Max_Iteration&&err>ERR);
	if(iteration==Max_Iteration)
		return 0;
	err=0;
	for(i=0;i<n[0];i++)
	{
		for(j=0;j<n[1];j++)
		{
			dh(i,j)*=sqrt(mu(i,j));
			h(i,j)+=dh(i,j);
			if(abs(dh(i,j))>err)
				err=abs(dh(i,j));
		}
	}
	printf("||dh||=%f\n",err);
	return iteration;
}

int set_image_data(BYTE* p,float *h,float *c,int offset,int type,int adjust)
{
	float temp,max=0,min=2000;;
	int i,j;
	if(type==1)             //save visible space
	{
		for(i=0;i<N[0];i++)
		{
			for (j=0;j<N[1];j++)
			{
				int i1,j1;
				temp=0;
				for(i1=0;i1<m[0];i1++)
					for(j1=0;j1<m[1];j1++)
						temp+=h(i+i1,j+j1)*c(m[0]-1-i1,m[1]-1-j1);
				if(temp>max)
					max=temp;
				if(temp<min)
					min=temp;
			}
		}
		if(max==min)
			max=min+1;
		for(i=0;i<N[0];i++)
		{
			for (j=0;j<N[1];j++)
			{
				int i1,j1;
				temp=0;
				for(i1=0;i1<m[0];i1++)
					for(j1=0;j1<m[1];j1++)
						temp+=h(i+i1,j+j1)*c(m[0]-1-i1,m[1]-1-j1);
				if(adjust)
					temp=(BYTE)(min+(temp-min)*(255-min)/(max-min));
				else
				{
					if(max>255)
						temp=(BYTE)(255*temp/max);
					else
						temp=(BYTE)(temp);
				}
				*(p+offset+3*(i*N[1]+j))=(BYTE)temp;
				*(p+offset+3*(i*N[1]+j)+1)=(BYTE)temp;
				*(p+offset+3*(i*N[1]+j)+2)=(BYTE)temp;
			}
		}
	}
	if(type==2)    //save hidden space
	{
		for(i=0;i<N[0];i++)
		{
			for (j=0;j<N[1];j++)
			{
				temp=h(i+(m[0]-1)/2,j+(m[1]-1)/2);
				if(temp>max)
					max=temp;
				if(temp<min)
					min=temp;
			}
		}
		if(max==min)
			max=min+1;
		for(i=0;i<N[0];i++)
			for (j=0;j<N[1];j++)
			{
				temp=h(i+(m[0]-1)/2,j+(m[1]-1)/2);
				if(adjust)
					temp=(BYTE)(min+(temp-min)*(255-min)/(max-min));
				else
				{
					if(max>255)
						temp=(255*temp/max);
					else
						temp=(temp);
				}
				*(p+offset+3*(i*N[1]+j))=(BYTE)temp;
				*(p+offset+3*(i*N[1]+j)+1)=(BYTE)temp;
				*(p+offset+3*(i*N[1]+j)+2)=(BYTE)temp;
			}
	}
	return 0;
}

int maximum_entropy(char *filename,char *p_directory,int mx,int my,float sigma_x,float sigma_y,BOOL test)
{
	char directory[200];
	if(create_directory_for_images(directory,p_directory,mx,my,psf.gamma,psf.pixel_to_distance)==0)
		return 0;

	int offset,size[2]={0,0};
	PVOID pimage=load_image(filename,size[0],size[1]);
	if(pimage==NULL)
	{
		printf("OpenFileError!\n");
		return 0;
	}

	N[0]=size[0];
	N[1]=size[1];        //as long as N[0]<=9000, memory will not overflow!
	m[0]=mx;
	m[1]=my;
	n[0]=m[0]+N[0]-1;
	n[1]=m[1]+N[1]-1;
	
	float *h,*h0,*dh,*mu,*g,*c,*b,*A,*B,*sigma;
	h=h0=dh=mu=g=c=b=A=B=sigma=NULL;
	if(alloc_initiate(h,h0,dh,mu,g,c,b,A,B,sigma)==0)
	{
		printf("Alloc memory error!\n");
		return 0;
	}
	
	//set_gray_pallette(pimage);
	offset=set_image_data_to_b(pimage,b,sigma);
	if(test)
		get_image(b,h,c,sigma);                       //this line is for testing, to construct original image
	//for(psf.gamma=0.8;psf.gamma<=1.2;psf.gamma+=0.1){
	//	for(omiga=0.12;omiga<=0.12;omiga+=0.05){
	set_Lorenze_probe(c);
	set_starting_point(h0);

	DWORD t1=GetTickCount();
	get_L_gradient2(A,c);
	DWORD t2=GetTickCount();
	printf("ellapsed time %d ms \n",t2-t1);

	float alfa=100,step,beta;
	int times=0,GS,GS_total=0;
	step=pow(1000,0.01);
	if(test)    //if testing, save the image automatically generated
	{
		set_image_data((BYTE*)pimage,h,c,offset,1,0);
		save_image_data(pimage,times,alfa,directory);
	}
	memcpy((PVOID)h,(PVOID)h0,sizeof(float)*n[0]*n[1]);

	do{
		times++;
		set_image_data((BYTE*)pimage,h,c,offset,2,1);
		save_image_data(pimage,times,alfa,directory);    //save images in every iteration, including starting image
		beta=alfa*2;
		memcpy((PVOID)mu,(PVOID)h,sizeof(float)*n[0]*n[1]);
		memcpy((PVOID)h0,(PVOID)h,sizeof(float)*n[0]*n[1]);
		get_g(g,h,h0,c,b,sigma,alfa);
		get_AB(A,B,c,mu,beta);
		GS=get_dh_SOR(B,mu,g,dh,h);
		GS_total+=GS;
		if(!GS)
		{
			printf("iteration cannot converge!\n");
			break;
		}
		printf("alfa=%f iteration=%d\n",alfa,GS);  
		alfa/=step;
	}while(alfa>0.1);

	printf("omiga=%f,alfa=%f, gamma=%f,iteration=%d\n",omiga,alfa,psf.gamma,GS_total);
	{
		float gamma=psf.gamma;
		psf.gamma=0.3;
		set_Lorenze_probe(c);
		set_image_data((BYTE*)pimage,h,c,offset,1,1);
		save_image_data(pimage,times+1,alfa-1,directory);
		psf.gamma=gamma;
	}
	//	}}
	if(alloc_free(h,h0,dh,mu,g,c,b,A,B,sigma)==0)
		printf("Free memory error!\n");
	VirtualFree(pimage,0,MEM_RELEASE);
	return 0;
};

int _tmain(int argc, _TCHAR* argv[])
{
	int i,j;
	float sigma=0,sigma_x=0,sigma_y=0;
	char filename[100]="";
	char directory[100]="";

	psf.Cs=1.2E7;
	psf.lamda=0.0197;
	psf.alfa=11./1000.;
	psf.df=-sqrt(psf.Cs*psf.lamda);
	psf.pixel_to_distance=3.581;


	printf("Enter file path and name:\n");
	//scanf("%s",filename);
	//if(strlen(filename)<=1)
		strcpy(filename,"F:\\MaxEnt\\test512.BMP");

	printf("Enter directory path to preserve images:\n");
	//scanf("%s",directory);
	//if(strlen(directory)<=1)
		strcpy(directory,"F:\\MaxEnt\\temp");

	printf("resolution:\n");
	//scanf_s("%f",&(psf.gamma));
	psf.gamma=2;
	psf.gamma/=2;

	printf("pixels/angstrom:\n");
	//scanf("%f",&(psf.pixel_to_distance));
	psf.pixel_to_distance=10;
	//float sita0,r0;
	//float chi,i_psf,i_gauss,chi_min=1E10,Cs_min,alfa_min,df_min,scale;
	//int dm=19,s,t;
	//sigma_x=5;
	//sigma_y=5;
	//for(psf.Cs=1.2E7;psf.Cs<=1.2E7;psf.Cs+=1E7)
	//{
	//	for(psf.alfa=11./1000.;psf.alfa<=11./1000.;psf.alfa+=1/1000.)
	//	{
	//		float df=-sqrt(psf.Cs*psf.lamda);
	//		for(psf.df=-df;psf.df>=df;psf.df-=10)
	//		{
	//			chi=0;
	//			for(i=(dm-1)/2;i<dm;i++)
	//			{
	//				for(j=(dm-1)/2;j<dm;j++)
	//				{
	//					s=i-(dm-1)/2;
	//					t=j-(dm-1)/2;
	//					sita0=0;
	//					r0=sqrt((float)(s*s+t*t))/psf.pixel_to_distance;
	//					i_psf=calculate_PSF(psf.Cs,psf.lamda,psf.alfa,psf.df,r0,sita0);
	//					i_gauss=exp(-(s*s/(sigma_x*sigma_x)+t*t/(sigma_y*sigma_y)));
	//					if(s==0&&t==0)
	//						scale=i_gauss/i_psf;
	//					chi+=(i_psf*scale-i_gauss)*(i_psf*scale-i_gauss);
	//				}
	//			}
	//			printf("Cs=%f,alfa=%f,df=%f,chi=%f\n",psf.Cs,psf.alfa,psf.df,chi);
	//			if(chi<chi_min)
	//			{
	//				chi_min=chi;
	//				Cs_min=psf.Cs;
	//				alfa_min=psf.alfa;
	//				df_min=psf.df;
	//			}
	//		}
	//	}
	//}
	//printf("Cs_min=%f,alfa_min=%f,df_min=%f,chi_min=%f\n",Cs_min,alfa_min,df_min,chi_min);
	for(i=19;i<=19;i+=2)
	{
		maximum_entropy(filename,directory,i,i,sigma_x,sigma_y,FALSE);
	}

	return 0;
}


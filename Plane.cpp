/**********************************************************************************************
*��������ˣ����ǳ�    ��ʼʱ�䣺2018.3.15��6.00    ���ʱ�䣺2018.3.17 20:45                 *
*�����ܣ�		(1).������֪������ƽ���һ�㷽��	                                          *
*				(2).����ƽ���һ�㷽�����ɷ������������Ϣ                                    *
*	       		(3).�ж�ĳ���Ƿ���ƽ����,��������õ���ƽ��ľ���[������Ĺ�ϵ]               *
*		   		(4).����ƽ���λ�ù�ϵ��������[������Ĺ�ϵ]                                *
*�����ŵ㣺 	(1).������ӱ��������Ŀ,����ǿ��                                               *
*				(2).�����˸������͵Ĺ��캯���ͳ�Ա����                                        *
*				(3).�������ô���������װ                                                      *
*				(4).д����ϸ��ע�ͣ������Ķ��͸Ľ�                                            *
*				(5).�漰�����ռ�������ε�֪ʶ,���������жϣ�system����...                  *
*				(6).�漰��ŷ��ɸ����������ɫ���ơ�������ʽ�����ȣ������Լ���ǰд�ģ���������  *
*����0(����)   (1).main()������     407-478��                                                 *
*              (2).Plane��Ķ��壺  80-181��                                                  *
*              (3).��Ա�����Ķ��壺 183-345��                                                 *
**********************************************************************************************/



#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <queue>
#include <string>
#include <vector>
#include <deque>
#include <cstdlib>
#include <stack>
#include <cctype>
#include <locale>
#include <iterator>
#include <sstream>
#include <windows.h>
#define N 3
#define a(i,j) a[i][j]
#define b(i,j) b[i][j]
#define RED    SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_RED)
#define BLUE   SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_BLUE)
#define GREEN  SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_GREEN)
#define WHITE  SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_RED|FOREGROUND_BLUE|FOREGROUND_GREEN)
#define CYAN   SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_BLUE|FOREGROUND_GREEN)
#define PURPLE SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_RED|FOREGROUND_BLUE)
#define YELLOW SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_RED|FOREGROUND_GREEN)
#define B_WHITE SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),BACKGROUND_INTENSITY|BACKGROUND_RED)
#define sl Sleep(100)
using namespace std;
#define For(i,a,b) for(int i=(a);i<=(b);++i)
#define Rep(i,a,b) for(int i=(a);i>=(b);--i)
#define _for(i,a,b) for(int i=(a);i<(b);++i)
#define _rep(i,a,b) for(int i=(a);i>(b);--i)
#define Mst(a,b) memset(a,(b),sizeof(a))
#define CL(a) memset(a,0,sizeof(a))
#define LL long long
#define ll long long
#define nullptr NULL
#define pb push_back
#define MP make_pair
#define pii pair<int,int>
#define PLL pair<long long,long long>
typedef struct xishu    //ƽ��һ�㷽�̵�ϵ��
{
	double A,B,C,D;
} Xs;
typedef	struct Vector   //�ռ������������
{
	double x,y,z,lenth;
	long long a,b;	//�ں��������ʱʹ��
} NomalVector;
typedef	struct coordinate   //�ռ������������
{
	double x,y,z,distance;
} Cd;
////////////////////////////////////  ����ȫ�ֱ��������㺯������  //
double point_3[3][3] = {0};                                      //
double a = 1.0,b = 1.0,c = 1.0,d =1.0;                           //
NomalVector n_vector,v1,v2;//����                                //
Cd point_coor,point_coor2,point_coor3;//������                              //
///////////////////////////////////////////////////////////////////
double Det(double a[N][N],int n);
double Dot_Product(double * a,double *b);
//�����ռ�������������

void   GetVectorProduct(double * a,double * b,double * vec);
//�����ռ�������������

void   GetVectorProduct(NomalVector v1,NomalVector v2,NomalVector & nv);

double Dot_Product(double x1=0,double y1=0,double z1=0,double x2=0,double y2=0,double z2=0);
//ͬ��

double norm_vector(double * a);
//��ռ�������ģ�����ɹ��캯��ʵ�ָù���

double norm_vector(double x=0,double y=0,double z=0);
//ͬ��

void   GetXishu(double a[3][3],double & s1,double &s2,double & s3,double &s0);

bool   IsPrime(int n);

////////////////////////////////////////////////////////////////
const int maxn = 1e7;
int xishu[700000];  //ָ��(ϵ��)
long long sushu[700000];
long long n[maxn];
void decomposition_prime_factor(long long number); //�ֽ�������
void Getsushubiao();
////////////////////////////////////////////////////////////////

long long gcd(long long a,long long b);
class Plane
{
	private:
		Xs p;   //ƽ��һ�㷽�̵�ϵ��
		NomalVector vec;  //����������
		NomalVector point; //ĳһ�������
		NomalVector point_on;
		void GetNormalVector();
		void GetPoint_on();
		bool IsOnThePlane(double x=1.0,double y=1.0,double z=1.0);
		void reduct();

	public:
		void SetPlane(double a=1,double b=1.0,double c=1.0,double d=1.0);
		void SetPlane(double dot[3][3]);       //��֪ƽ���ϲ����ߵ�����ȷ�������÷���ϵ��
		void SetPlane(NomalVector nv,Cd cd1);       //��֪ƽ����һ��ͷ�����ȷ�������÷���ϵ��
		void SetPlane(NomalVector v1,NomalVector v2,Cd point);   //��֪ƽ����һ���������������ȷ�������÷���ϵ��
		void SetPlane(NomalVector v1,Cd p1,Cd p2);        //��֪ƽ���������һ����������ȷ�������÷���ϵ��
		void ShowAllInformation();
		Plane(double a=1,double b=1.0,double c=1.0,double d=1.0) //���캯��
		{
			if(a||b||c)
			{
				this->p.A = a;
				this->p.B = b;
				p.C = c;
				p.D = d;
				reduct();
			}
			GetNormalVector();
			GetPoint_on();

		}
		Plane()
		{
			//�޲ι��캯��
			this->p.A = 1.0;
			this->p.B = 1.0;
			this->p.C = 1.0;
			this->p.D = 0.0;
			GetNormalVector();
			GetPoint_on();
			//reduct();        //�޲ι��죬����Լ�֣�Ԥ��Ķ���1.0
		}
		Plane(const Plane &Q)
		{
			//�������캯��
			this->p.A = Q.p.A;
			p.B = Q.p.B;
			p.C = Q.p.C;
			p.D = Q.p.D;
			GetNormalVector();
			GetPoint_on();
			reduct();
		}
		Plane(double dot[3][3])                  //���캯������֪ƽ���ϲ����ߵ�����ȷ������ϵ��
		{
			GetXishu(dot,p.A,p.B,p.C,p.D);
			if(p.A==0&&p.B==0&&p.C==0)    //��֤��ƽ�淽�̣�A��B��C����ͬʱΪ��
			{
				p.A = p.B = p.C =1.0;
			}
			reduct();
			GetNormalVector();
			GetPoint_on();

		}
		Plane(NomalVector nv,Cd cd1)
		{
			//���캯������֪ƽ����һ��ͷ�����ȷ������ϵ��
			if(nv.x||nv.y||nv.z)
			{
				this->p.A = nv.x;
				this->p.B = nv.y;
				this->p.C = nv.z;
				this->p.D = -(nv.x*cd1.x+nv.y*cd1.y+nv.z*cd1.z);
			}
			else
			{
				this->p.A = p.B = p.C = p.D = 1.0;
			}
			reduct();
			GetNormalVector();
			GetPoint_on();
			ShowPlane();
		}
		Plane(NomalVector v1,NomalVector v2,Cd point)
		{
			//���캯������֪ƽ����һ���������������ȷ������ϵ��
			NomalVector nvector;
			GetVectorProduct(v1,v2,nvector);
			if(nvector.x||nvector.y||nvector.z)
			{
				this->p.A = nvector.x;
				this->p.B = nvector.y;
				this->p.C = nvector.z;
				this->p.D = -(nvector.x*point.x + nvector.y*point.y + nvector.z*point.z);
			}
			else
			{
				this->p.A = p.B = p.C = p.D = 1.0;
			}
			reduct();
			GetNormalVector();
			GetPoint_on();
			ShowPlane();

		}
		Plane(NomalVector v1,Cd p1,Cd p2)               //���캯������֪ƽ���������һ����������ȷ������ϵ��
		{
			NomalVector nvector,v2;
			v2.x = p1.x-p2.x;
			v2.y = p1.y-p2.y;
			v2.z = p1.z-p2.z;
			GetVectorProduct(v1,v2,nvector);
			if(nvector.x||nvector.y||nvector.z)
			{
				this->p.A = nvector.x;
				this->p.B = nvector.y;
				this->p.C = nvector.z;
				this->p.D = -(nvector.x*p1.x + nvector.y*p1.y + nvector.z*p1.z);
			}
			else
			{
				this->p.A = p.B = p.C = p.D = 1.0;
			}
			GetNormalVector();
			GetPoint_on();
			reduct();
		}
		void ShowPlane();
		void ShowNormalVector();
		void ShowPoint_on();
		void GetNormalVectorLenth();
		void JudgePointOnThePlane();
		friend void GetXishu(double a[3][3],double & s1,double &s2,double & s3,double &s0);
		friend void GetVectorProduct(NomalVector v1,NomalVector v2,NomalVector & nv);
		friend long long gcd(long long a,long long b);
		friend void pos_relation_between_plane(Plane p1,Plane p2);
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Plane::SetPlane(double a,double b,double c,double d)
{
	if(a||b||c)
	{
		this->p.A = a;
		this->p.B = b;
		p.C = c;
		p.D = d;
		reduct();
		GetNormalVector();
		GetPoint_on();
	}
	else
	{
		cout<<"���ݲ��Ϸ���A,B,C����ͬʱΪ�㣡\n";
	}
}
void Plane::SetPlane(double dot[3][3])      //��֪ƽ���ϲ����ߵ�����ȷ�������÷���ϵ��
{
	GetXishu(dot,p.A,p.B,p.C,p.D);
	if(p.A==0&&p.B==0&&p.C==0)    //��֤��ƽ�淽�̣�A��B��C����ͬʱΪ��
	{
		p.A = p.B = p.C =1.0;
	}
	reduct();
	GetNormalVector();
	GetPoint_on();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Plane::SetPlane(NomalVector nv,Cd cd1)           //��֪ƽ����һ��ͷ�����ȷ�������÷���ϵ��
{
	if(nv.x||nv.y||nv.z)
	{
		this->p.A = nv.x;
		this->p.B = nv.y;
		this->p.C = nv.z;
		this->p.D = -(nv.x*cd1.x+nv.y*cd1.y+nv.z*cd1.z);
	}
	else
	{
		this->p.A = p.B = p.C = p.D = 1.0;
	}
	reduct();
	GetNormalVector();
	GetPoint_on();
	ShowPlane();
}
void Plane::SetPlane(NomalVector v1,NomalVector v2,Cd point)        //��֪ƽ����һ���������������ȷ�������÷���ϵ��
{
	NomalVector nvector;
	GetVectorProduct(v1,v2,nvector);
	if(nvector.x||nvector.y||nvector.z)
	{
		this->p.A = nvector.x;
		this->p.B = nvector.y;
		this->p.C = nvector.z;
		this->p.D = -(nvector.x*point.x + nvector.y*point.y + nvector.z*point.z);
	}
	else
	{
		this->p.A = p.B = p.C = p.D = 1.0;
	}
	reduct();
	GetNormalVector();
	GetPoint_on();
	ShowPlane();
}
void Plane::SetPlane(NomalVector v1,Cd p1,Cd p2)            //��֪ƽ���������һ����������ȷ�������÷���ϵ��
{
	NomalVector nvector,v2;
	v2.x = p1.x-p2.x;
	v2.y = p1.y-p2.y;
	v2.z = p1.z-p2.z;
	GetVectorProduct(v1,v2,nvector);
	if(nvector.x||nvector.y||nvector.z)
	{
		this->p.A = nvector.x;
		this->p.B = nvector.y;
		this->p.C = nvector.z;
		this->p.D = -(nvector.x*p1.x + nvector.y*p1.y + nvector.z*p1.z);
	}
	else
	{
		this->p.A = p.B = p.C = p.D = 1.0;
	}
	reduct();
	GetNormalVector();
	GetPoint_on();

}
void Plane::ShowPlane()
{
	cout<<"��ƽ��ķ���Ϊ��\n";
	cout<<"("<<p.A<<")X + ("<<p.B<<")Y + ("<<p.C<<")Z + ("<<p.D<<") = 0\n";
	GREEN;
	if(p.A)  //A!=0
	{
		if(p.A==1.0) cout<<"X";
		else cout<<p.A<<"X";
		if(p.B) //A:=0,B!=0
		{
			if(p.B==1.0) cout<<" + Y";
			else if(p.B==-1.0) cout<<" - Y";
			else if(p.B>0) cout<<" + "<<p.B<<"Y";
			else cout<<" - "<<-p.B<<"Y";
		}
	}
	else  //A:=0
	{
		if(p.B)
		{
			if(p.B==1.0) cout<<"Y";
			else if(p.B==-1.0) cout<<"-Y";
			else cout<<p.B<<"Y";
		}
	}
	if(p.A||p.B)
	{
		if(p.C)
		{
			if(p.C==1.0) cout<<" + "<<"Z";
			else if(p.C==-1.0) cout<<" - "<<"Z";
			else if(p.C>0) cout<<" + "<<p.C<<"Z";
			else cout<<" - "<<-p.C<<"Z";
		}
	}
	else
	{
		if(p.C==1) cout<<"Z";
		else if(p.C==-1) cout<<"-Z";
		else cout<<p.C<<"Z";
	}
	if(p.D)
	{
		if(p.D>0) cout<<" + "<<p.D;
		else cout<<" - "<<-p.D;
	}
	cout<<" = 0\n";
	WHITE;
}
void Plane::GetNormalVector()
{
	vec.x = p.A;
	vec.y = p.B;
	vec.z = p.C;
	vec.lenth = sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
}
void Plane::GetPoint_on()
{
	double a = p.A,b = p.B,c = p.C,d = p.D;
	if(c==0&&a!=0) swap(a,c);
	if(c==0&&a==0) swap(b,c);
	point_on.x = 0.0;
	point_on.y = 1.0;
	point_on.z = -(b+d)/c;
	//cout<<"point_on.a = "<<-(p.B+p.D)*10000<<"  point_on.b = "<<p.C*10000<<endl;
	point_on.a = -(b+d)*10000;
	point_on.b = c*10000;
	long long max_common_div = gcd(point_on.a,point_on.b);
	//cout<<"max_common_div = "<<max_common_div<<endl;
	point_on.a /= max_common_div;
	point_on.b /= max_common_div;
	if(point_on.a>0&&point_on.b<0)
	{
		point_on.a *= -1;
		point_on.b *= -1;
	}
	//cout<<"point_on.a = "<<point.a<<" point_on.b = "<<point_on.b<<endl;
}
void Plane::ShowPoint_on()
{
	cout<<"ƽ����ĳһ�������Ϊ��("<<point_on.x<<" , "<<point_on.y<<" , ";
	printf("%.8lf )\n",point_on.z);//<<setprecision(10)<<point_on.z<<")\n";
	cout<<"�÷�����ʾ��һ��Ϊ��  ("<<point_on.x<<" , "<<point_on.y<<" , "<<point_on.a;
	if(point_on.b!=1) cout<<"/"<<point_on.b;
	cout<<" )\n";
}
void Plane::ShowNormalVector()
{
	cout<<"��ƽ��ķ���������Ϊ��("<<vec.x<<","<<vec.y<<","<<vec.z<<")\n";
}
bool Plane::IsOnThePlane(double x,double y,double z)
{
	double ss = fabs(x*p.A + y*p.B + z*p.C + p.D);
	point.lenth = ss/(sqrt(p.A*p.A+p.B*p.B+p.C*p.C));
	if(ss<1e-6) return true;   //����
	else return false;
}
void Plane::ShowAllInformation()
{
	cout<<"----------------------------------------------------------------------\n��ƽ����Ϣ���£�\n";
	ShowPlane();
	ShowNormalVector();
	ShowPoint_on();
	cout<<"----------------------------------------------------------------------\n";
}
void Plane::JudgePointOnThePlane()
{
	WHITE;
	cout<<"�����ж�ĳһ���Ƿ��ڸ�ƽ���ϣ�������������(x,y,z):\n";
	CYAN;
	cin>>point.x>>point.y>>point.z;
	WHITE;
	if(IsOnThePlane(point.x,point.y,point.z))
		cout<<"��("<<point.x<<","<<point.y<<","<<point.z<<")�ڸ�ƽ���ϡ�\n";
	else
	{
		cout<<"��("<<point.x<<","<<point.y<<","<<point.z<<")���ڸ�ƽ���ϡ�����ƽ��ľ���Ϊ��"<<point.lenth;
		cout<<"("<<fabs(point.x*p.A+point.y*p.B+point.z*p.C+p.D)<<"/sqr("<<p.A*p.A+p.B*p.B+p.C*p.C<<"))\n";
		GREEN;
		decomposition_prime_factor((long long)(p.A*p.A+p.B*p.B+p.C*p.C) );
		WHITE;
	}
}
void pos_relation_between_plane(Plane p1,Plane p2)
{
	WHITE;
	cout<<"������д��ƽ��ķ��̣�\n";
	YELLOW;
	cout<<"##################################################\n";
	PURPLE;
	cout<<"#Plane1��";
	p1.ShowPlane();
	PURPLE;
	cout<<"#Plane2��";
	p2.ShowPlane();
	YELLOW;
	cout<<"##################################################\n";
	WHITE;
	int tap = 0,i,j;
	bool flag = true;
	double a[5]= {0},b[5]= {0}; //double record[4]={0};
	a[1] = p1.p.A;
	a[2] = p1.p.B;
	a[3] = p1.p.C;
	a[4] = p1.p.D;
	b[1] = p2.p.A;
	b[2] = p2.p.B;
	b[3] = p2.p.C;
	b[4] = p2.p.D;
	if( a[1]*b[2]!=b[1]*a[2] || a[1]*b[3]!=b[1]*a[3] || a[2]*b[3]!=b[2]*a[3] ) flag = false;
	for(i=1; i<4; i++)
	{
		if(a[i]&&b[i]==0.0||b[i]&&a[i]==0.0)
		{
			flag = false;
			break;
		}
	}
	if(!flag) cout<<"��ƽ�治ƽ�С�\n";
	if(flag)
	{
		if(a[1]*b[2]==b[1]*a[2] && a[1]*b[3]==b[1]*a[3] && a[2]*b[3]==b[2]*a[3] &&
		        a[4]*b[1]==b[4]*a[1] && a[4]*b[2]==b[4]*a[2] && a[4]*b[3]==b[4]*a[3] )
			cout<<"��ƽ���غϣ�\n";
		else cout<<"��ƽ��ƽ�С�\n";
	}
	if(flag) return;
	double x[3],y[3];
	x[0] = p1.vec.x;
	x[1] = p1.vec.y;
	x[2] = p1.vec.z;
	y[0] = p2.vec.x;
	y[1] = p2.vec.y;
	y[2] = p2.vec.z;
	double pot_produ = fabs(Dot_Product(x,y));
	double mol_produ = (x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) * (y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
	double cos_angle = pot_produ/(sqrt(mol_produ));
	double angle = acos(cos_angle);
	cout<<"��ƽ���������ǵ�����ֵΪ"<<setprecision(8)<<cos_angle<<" ("<<pot_produ<<"/(sqr("<<mol_produ<<"))\n";
	GREEN;
	decomposition_prime_factor((long long) (mol_produ));
	WHITE;
	cout<<"����ǵĻ���ֵԼΪ��"<<angle<<"(ԼΪ"<<angle*180<<"��)\n";
	cout<<"��ƽ�潻����������ƽ�淽�������ɵķ����顣\n";


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void showfunction(void);
void function_1(Plane & p);
void function_2(Plane & p);
void function_3(Plane & p);
void function_4(Plane & p);
void function_5(Plane & p);
void function_6(Plane & p1,Plane & p2);
void fun_point(Plane p);
int main()
{
	system("title ƽ�����by���ǳ�");
	int i,j,choose;
	char flag = '1';
	CYAN;
	cout<<"���ǳ�������2018.3.15-3.17\n";
	WHITE;
	cout<<"ƽ���һ�㷽������ Ax + By + Cz + D = 0;\n";
	PURPLE;
	cout<<"�����������¹��ܣ�\n1.������֪������ƽ���һ�㷽��\n2.����ƽ���һ�㷽�����ɷ������������Ϣ\n";
	cout<<"3.�ж�ĳ���Ƿ���ƽ����,��������õ���ƽ��ľ���\n";
	Getsushubiao();
	Plane p1(0,-1,2,-5);
	Plane p2(p1);

	int flagg = 0; //��֤��ִ��һ�β�������
	while(flag!='0')
	{
		if(flag!='1'&&flag!='2') goto go_on;
		if(flagg&&flag=='2')
		{
			system("CLS"); //����
			CYAN;
			cout<<"���ǳ�������\n";
		}
		showfunction();
		cout<<"������1~6�����֣�_\b";
Line:
		cout<<"_\b";
		CYAN;
		cin>>choose;
		WHITE;
		switch(choose)
		{
			case 1:
				function_1(p1);
				break;
			case 2:
				function_2(p1);
				break;
			case 3:
				function_3(p1);
				break;
			case 4:
				function_4(p1);
				break;
			case 5:
				function_5(p1);
				break;
			case 6:
				function_6(p1,p2);
				break;
defalt:
				cout<<"����������������룡\n";
				goto Line;
		}
go_on:
		YELLOW;
		cout<<"����ʹ���밴";
		GREEN;
		cout<<"1";
		YELLOW;
		cout<<"����ʹ�ò������밴";
		PURPLE;
		cout<<"2";
		YELLOW;
		cout<<"�˳��밴";
		RED;
		cout<<"0";
		YELLOW;
		cout<<"��_\b";
		CYAN;
		cin>>flag;
		flagg = 1;
		cout<<endl<<endl;
	}
	cout<<"\n\n";
	char ch[] = "Bye-bye!See you next time!";
	for(i=0; ch[i]!='\0'; i++)
	{
		switch(i%7)
		{
			case 0:
				CYAN;
				break;
			case 1:
				BLUE;
				break;
			case 2:
				PURPLE;
				break;
			case 3:
				GREEN;
				break;
			case 4:
				YELLOW;
				break;
			case 5:
				RED;
				break;
			case 6:
				WHITE;
				break;
		}
		cout<<ch[i];
		sl;
	}
	WHITE;
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Det(double a[N][N],int n)//����ʽ��������)
{
	if(n==1)return a[0][0];
	if(n==2)return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
	if(n>2)
	{
		double b[N][N]= {0},answer=0.0;
		int i,j,k,l;
		for(i=0; i<n; i++)
		{
			for(j=0; j<n-1; j++)
			{
				for(k=0; k<n-1; k++)
				{
					if(k<i) b[j][k]=a[j+1][k];
					else    b[j][k]=a[j+1][k+1];

				}
			}
			answer+=Det(b,n-1)*pow(-1,i)*a[0][i];
		}
		return answer;
	}
}

double Dot_Product(double * a,double *b)
{
	double s = 0.0;
	for(int i=0; i<3; i++)
		s += a[i]*b[i];
	return s;
}
double Dot_Product(double x1,double y1,double z1,double x2,double y2,double z2)
{
	return x1*x2 + y1*y2 + z1*z2;
}
double norm_vector(double * a)
{
	return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
double norm_vector(double x,double y,double z)
{
	return sqrt(x*x + y*y + z*z);
}
void GetXishu(double a[3][3],double & s1,double &s2,double & s3,double &s0)  //Plane�����Ԫ����
{
	double b(3,3) = {0};
	int i,j;
	s0 = Det(a,3);
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
			b(i,j) = a(i,j);
	}
	b(0,0) = b(1,0) = b(2,0) = 1.0;
	s1 = -Det(b,3);
	b(0,1) = a(0,0);
	b(1,1) = a(1,0);
	b(2,1) = a(2,0);
	s2 = Det(b,3);
	///////////////////////////////////////////////////////////
	b(0,2) = a(0,1);
	b(1,2) = a(1,1);
	b(2,2) = a(2,1);
	//////////////////////////////////
	s3 = -Det(b,3);
}
void GetVectorProduct(double * a,double * b,double * vec)
{
	vec[0] = a[1]*b[2] - a[2]*b[1];
	vec[1] = b[0]*a[2] - a[0]*b[2];
	vec[2] = b[1]*a[0] - a[1]*b[0];
}
void GetVectorProduct(NomalVector v1,NomalVector v2,NomalVector & nv)
{
	nv.x = v1.y*v2.z - v1.z*v2.y;
	nv.y = v1.z*v2.x - v1.x*v2.z;
	nv.z = v1.x*v2.y - v1.y*v2.x;
}

int gcd(int a,int b)
{
	return b==0?a:gcd(b,a%b);
}
bool IsPrime(int n)
{
	if(n==1||n==0) return false;
	for(int i=2; i<n; i++)
	{
		if(n%i==0) return false;
	}
	return true;
}
long long gcd(long long a,long long b)
{
	return b==0?a:gcd(b,a%b);
}
void Plane::reduct()  //���̻���ͬʱ����1000����Լ��
{
	bool F[4] = {1,1,1,1};
	int i,j,num[4]= {0},minn,gcdmax(1);
	long long a[4] = {0};
	a[0] = this->p.A*1000;
	a[1] = this->p.B*1000;
	a[2] = this->p.C*1000;
	a[3] = this->p.D*1000;
	for(i=0; i<4; i++)
	{
		if(a[i]) minn = abs(a[i]);
	}
	for(i=0; i<4; i++)
	{
		num[i] = a[i];
		if(a[i]<0)
		{
			num[i]*=-1;
			F[i] = false;
		}
		if(minn>num[i]&&num[i]) minn = num[i];
	}
	//cout<<"minn = "<<minn<<endl;
	for(i=0; i<4; i++)
	{
		if(IsPrime(abs(a[i])))
		{
			goto Line;
		}
	}
	for(i=minn; i>0; i--)
	{
		for(j=0; j<4; j++)
		{
			if(num[j]%i) break;
		}
		if(j==4)
		{
			gcdmax = i;
			break;
		}
	}
Line:
	//cout<<num[0]<<" "<<num[1]<<" "<<num[2]<<" "<<num[3];
	//cout<<"�⼸���������Լ���ǣ�"<<gcdmax<<endl;
	if(!F[0]) gcdmax *=-1;
	for(i=0; i<4; i++)
	{
		num[i] /= gcdmax;
		if(!F[i]) num[i] *=-1;
	}

	this->p.A = (double)num[0];
	this->p.B = (double)num[1];
	this->p.C = (double)num[2];
	this->p.D = (double)num[3];
	if(p.A==0&&p.B<0)                    //ʹ���̵ĵ�һ��������ϵ��Ϊ��
	{
		p.B *= -1;
		p.C *= -1;
		p.D *= -1;
	}
	if(p.A==0&&p.B==0&&p.C<0)            //ͬ��
	{
		p.C *= -1;
		p.D *= -1;
	}

}
void showfunction(void)
{
	YELLOW;
	cout<<"-----------------------------------------------------------------------------------------------------\n";
	cout<<"|(1).��֪ƽ��һ�㷽�̵�ϵ���밴1                 (2).��֪ƽ���ϲ����ߵ�������������ƽ�淽���밴2    \n";
	cout<<"|(3).��֪ƽ����һ���������ƽ�淨�����󷽳��밴3 (4).��֪ƽ����һ���������������������󷽳��밴4    \n";
	cout<<"|(5).��֪ƽ�������������һ�����������󷽳��밴5	(6).�ж�ƽ���λ�ù�ϵ���������밴6\n";
	cout<<"-----------------------------------------------------------------------------------------------------\n";
	WHITE;
}

void function_1(Plane & p)
{
	cout<<"�����뷽�̵�ϵ����A,B,C,D��_ _ _ _\b\b\b\b\b\b\b";
	CYAN;
	cin>>a>>b>>c>>d;
	WHITE;
	p.SetPlane(a,b,c,d);
	p.ShowAllInformation();
	fun_point(p);
}

void function_2(Plane & p)
{
	int i,j;
	cout<<"���ǿ���ͨ��ĳƽ���ϲ����ߵ�������ȷ�����ƽ�档����������������꣺\n";
	for(i=0; i<3; i++)
	{
		WHITE;
		cout<<"_ _ _\b\b\b\b\b";
		CYAN;
		for(j=0; j<3; j++)
			cin>>point_3[i][j];
	}
	WHITE;
	p.SetPlane(point_3);
	p.ShowAllInformation();
	fun_point(p);
}
void function_3(Plane & p)
{
	cout<<"���ǻ�����ͨ��ƽ����һ�������ƽ��ķ�����ȷ��ƽ�棺\n�����뷨�������꣺_ _ _\b\b\b\b\b";
	CYAN;
	cin>>n_vector.x>>n_vector.y>>n_vector.z;
	WHITE;
	cout<<"������ƽ��ĳ������꣺_ _ _\b\b\b\b\b";
	CYAN;
	cin>>point_coor.x>>point_coor.y>>point_coor.z;
	WHITE;
	p.SetPlane(n_vector,point_coor);
	p.ShowAllInformation();
	fun_point(p);
}
void function_4(Plane & p)
{
	cout<<"������ƽ������֪�ĵ����꣺_ _ _\b\b\b\b\b";
	CYAN;
	cin>>point_coor.x>>point_coor.y>>point_coor.z;
	WHITE;
	cout<<"�������һ���������������꣺_ _ _\b\b\b\b\b";
	CYAN;
	cin>>v1.x>>v1.y>>v1.z;
	WHITE;
	cout<<"������ڶ����������������꣺_ _ _\b\b\b\b\b";
	CYAN;
	cin>>v2.x>>v2.y>>v2.z;
	WHITE;
	p.SetPlane(v1,v2,point_coor);
	p.ShowAllInformation();
	fun_point(p);
}
void function_5(Plane & p)
{
	cout<<"���ǻ�����ͨ����֪ƽ���������һ����������ȷ������ϵ����\n";
	cout<<"������ƽ��ķ������������꣺_ _ _\b\b\b\b\b";
	CYAN;
	cin>>v1.x>>v1.y>>v1.z;
	WHITE;
	cout<<"������ƽ������֪�ĵ�һ�������꣺_ _ _\b\b\b\b\b";
	CYAN;
	cin>>point_coor.x>>point_coor.y>>point_coor.z;
	WHITE;
	cout<<"������ƽ������֪�ĵڶ��������꣺_ _ _\b\b\b\b\b";
	CYAN;
	cin>>point_coor2.x>>point_coor2.y>>point_coor2.z;
	WHITE;
	p.SetPlane(v1,point_coor,point_coor2);
	p.ShowAllInformation();
	fun_point(p);
}

void function_6(Plane & p1,Plane & p2)
{
	int n;
	PURPLE;
	cout<<"=====================================================================================================\n";
	WHITE;
	cout<<"�ж���ƽ���λ�ù�ϵ������Ҫ��ȡ��ƽ��ķ��̡�\n";
	showfunction();
	cout<<"������1~5����������ȡ��һ��ƽ�����Ϣ(��ֻ�õ����̼���)��_\b";
	CYAN;
	cin>>n;
	WHITE;
	switch(n)
	{
		case 1:
			function_1(p1);
			break;
		case 2:
			function_2(p1);
			break;
		case 3:
			function_3(p1);
			break;
		case 4:
			function_4(p1);
			break;
		case 5:
			function_5(p1);
			break;
	}
	WHITE;
	cout<<"����ͼ��������1~5����������ȡ��һ��ƽ�����Ϣ(��ֻ�õ����̼���)��_\b";
	CYAN;
	cin>>n;
	WHITE;
	switch(n)
	{
		case 1:
			function_1(p2);
			break;
		case 2:
			function_2(p2);
			break;
		case 3:
			function_3(p2);
			break;
		case 4:
			function_4(p2);
			break;
		case 5:
			function_5(p2);
			break;
	}
	pos_relation_between_plane(p1,p2);
	PURPLE;
	cout<<"=====================================================================================================\n";
	WHITE;
}

void fun_point(Plane p)
{
	char flag = '1';
	WHITE;
	cout<<"��ĳ�㵽ƽ��ľ����밴1�������밴0��_\b";
	CYAN;
	cin>>flag;
	while(flag!='0')
	{
		if(flag!='1')  goto Line1;
		p.JudgePointOnThePlane();
Line1:
		WHITE;
		cout<<"�����������������ƽ��ľ����밴1�������밴0��_\b";
		CYAN;
		cin>>flag;
		WHITE;
	}
}


void Getsushubiao()
{
	int i,j;
	for(i=0; i<maxn; i++)   //��׼������������ֵ���±���ͬ
	{
		n[i] = i;
	}
	for(i=2; i<maxn; i++)  //ŷ��ɸ������
	{
		for(j=i*2; j<maxn; j+=i)
		{
			n[j] = 0;
		}
	}
	j = -1;
	for(i=2; i<maxn; i++) //��ɸ�õ���������������ֽ�������
	{
		if(n[i])
			sushu[++j] = n[i];
	}
}

void decomposition_prime_factor(long long number)
{
	int i,j;
	long long num = number;
	if(!number) return;
	memset(xishu,0,sizeof(xishu));
	for(i=0; i<664578; i++)
	{
		if(num==1) break;
		if(num<sushu[i]) break;
		while(num%sushu[i]==0)
		{
			num/=sushu[i];
			xishu[i]++;
		}
	}
	printf("(PS:�ֽ�������,����л��^_^)  %lld = ",number);
	for(int k=0; k<i-1; k++)
	{
		if(!xishu[k]) continue;
		if(xishu[k]>1) printf("%d^%d * ",sushu[k],xishu[k]);
		else printf("%d * ",sushu[k]);
	}
	if(number==1) printf("2^0\n");
	else if(xishu[i-1]>1)
		printf("%d^%d\n",sushu[i-1],xishu[i-1]);
	else printf("%d\n",sushu[i-1]);
}

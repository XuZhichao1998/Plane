#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cstring>
using namespace std;
int main(void)
{
	double a[10]={0},b[10]={0},vector_product[10]={0};//������ 
	cout<<"����ռ�����������������������\n";
	int flag(1),i(0);
	double dot_product = 0.0;//��ˣ��������� 
	while(flag)
	{
		memset(a,0,sizeof(a));
		memset(b,0,sizeof(b));
		memset(vector_product,0,sizeof(vector_product));
		dot_product = 0.0;
		printf("����������a�ĺ����ꡢ�����ꡢ�����꣺\n");
		for(i=0;i<3;i++)
			scanf("%lf",a+i);
		printf("����������b�ĺ����ꡢ�����ꡢ�����꣺\n");
		for(i=0;i<3;i++)
			scanf("%lf",b+i);
		cout<<"����a:"<<"("<<a[0]<<","<<a[1]<<","<<a[2]<<")"<<endl;
		cout<<"����b:"<<"("<<b[0]<<","<<b[1]<<","<<b[2]<<")"<<endl;
		for(i=0;i<3;i++)
			dot_product += a[i] * b[i]; 	
		vector_product[0] = a[1]*b[2] - a[2]*b[1];
		vector_product[1] = b[0]*a[2] - a[0]*b[2];
		vector_product[2] = b[1]*a[0] - a[1]*b[0];
		cout<<"����a��b��������Ϊ��"<<dot_product<<endl;
		cout<<"����a��b������������Ϊ��("<<vector_product[0]<<","<<vector_product[1]<<","<<vector_product[2]<<")\n";
		cout<<"\n����ʹ��������1���˳�������0��";
		cin>>flag; 	
	}  
	
	return 0;
}

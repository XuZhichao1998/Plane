#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cstring>
using namespace std;
int main(void)
{
	double a[10]={0},b[10]={0},vector_product[10]={0};//向量积 
	cout<<"计算空间向量的数量积、向量积：\n";
	int flag(1),i(0);
	double dot_product = 0.0;//点乘（数量积） 
	while(flag)
	{
		memset(a,0,sizeof(a));
		memset(b,0,sizeof(b));
		memset(vector_product,0,sizeof(vector_product));
		dot_product = 0.0;
		printf("请输入向量a的横坐标、纵坐标、竖坐标：\n");
		for(i=0;i<3;i++)
			scanf("%lf",a+i);
		printf("请输入向量b的横坐标、纵坐标、竖坐标：\n");
		for(i=0;i<3;i++)
			scanf("%lf",b+i);
		cout<<"向量a:"<<"("<<a[0]<<","<<a[1]<<","<<a[2]<<")"<<endl;
		cout<<"向量b:"<<"("<<b[0]<<","<<b[1]<<","<<b[2]<<")"<<endl;
		for(i=0;i<3;i++)
			dot_product += a[i] * b[i]; 	
		vector_product[0] = a[1]*b[2] - a[2]*b[1];
		vector_product[1] = b[0]*a[2] - a[0]*b[2];
		vector_product[2] = b[1]*a[0] - a[1]*b[0];
		cout<<"向量a和b的数量积为："<<dot_product<<endl;
		cout<<"向量a和b的向量积坐标为：("<<vector_product[0]<<","<<vector_product[1]<<","<<vector_product[2]<<")\n";
		cout<<"\n继续使用请输入1，退出请输入0：";
		cin>>flag; 	
	}  
	
	return 0;
}

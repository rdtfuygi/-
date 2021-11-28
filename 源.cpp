#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
using namespace std;

#define Ants 100			//蚂蚁的个数
#define C 21				//地点的数量
#define I 50				//最大迭代次数
#define Alpha 1				//表征信息素重要程度的参数
#define Beta 5				//表征启发式因子重要程度的参数
#define Rho 0.1				//信息素蒸发系数
#define Q 100				//信息素增加强度系数

string qi_dian_b[2] = { "南食","北食" };
string di_dian[C] = { "南食","南1","南2","南3","南4","南5","南6","南7","南8","南9","南10" "北1","北2","北3","北4","北5","北6","北7","北8","北9","北10","北11" };
int lu_jing_dian_shu = 0;
int lu_jing_dian[C];


int Tabu[Ants][C];			//禁忌表，存储走过的路径
int Dist[C][C] =		//表示两两地点间距离,直接用步行时间表示
{
	{0,13,12,1,6,6,2,14,5,6,4,1,7,2,12,7,10,2,3,3,15},
	{13,0,7,1,5,15,13,2,8,3,2,3,12,12,14,12,2,1,12,14,15},
	{12,7,0,4,8,12,7,11,14,8,11,6,5,11,9,3,14,5,15,4,14},
	{1,1,4,0,2,3,11,2,15,5,15,10,13,4,13,11,8,2,7,3,4},
	{6,5,8,2,0,4,5,3,3,9,12,13,8,15,2,15,14,3,12,8,15},
	{6,15,12,3,4,0,15,5,15,11,15,13,10,13,3,12,7,1,5,8,3},
	{2,13,7,11,5,15,0,15,7,7,13,7,7,6,8,13,3,5,6,10,9},
	{14,2,11,2,3,5,15,0,8,13,5,7,10,1,3,1,9,14,13,3,7},
	{5,8,14,15,3,15,7,8,0,14,11,3,2,8,9,5,6,15,14,14,13},
	{6,3,8,5,9,11,7,13,14,0,12,3,12,7,2,14,1,8,9,2,6},
	{4,2,11,15,12,15,13,5,11,12,0,3,2,10,6,13,13,11,5,3,9},
	{1,3,6,10,13,13,7,7,3,3,3,0,11,13,8,9,2,1,6,10,9},
	{7,12,5,13,8,10,7,10,2,12,2,11,0,13,13,2,12,13,14,11,15},
	{2,12,11,4,15,13,6,1,8,7,10,13,13,0,15,15,12,11,9,4,14},
	{12,14,9,13,2,3,8,3,9,2,6,8,13,15,0,8,11,9,11,9,12},
	{7,12,3,11,15,12,13,1,5,14,13,9,2,15,8,0,11,5,11,14,13},
	{10,2,14,8,14,7,3,9,6,1,13,2,12,12,11,11,0,2,12,10,7},
	{2,1,5,2,3,1,5,14,15,8,11,1,13,11,9,5,2,0,9,11,5},
	{3,12,15,7,12,5,6,13,14,9,5,6,14,9,11,11,12,9,0,6,9},
	{3,14,4,3,8,8,10,3,14,2,3,10,11,4,9,14,10,11,6,0,2},
	{15,15,14,4,15,3,9,7,13,6,9,9,15,14,12,13,7,5,9,2,0}
};
double Eta[C][C];			//表示启发式因子，为Dist[][]中距离的倒数
double Tau[C][C];			//表示两两节点间信息素浓度
double DeltaTau[C][C];		//表示两两节点间信息素的变化量
double L_best[I];			//存储每次迭代的路径的最短长度
double L_ave[I];			//存储每次迭代的路径的平均长度
int R_best[I][C];			//存储每次迭代的最佳路线

void hb()
{
	int d2[C][C] =
	{
		{0,13,10,14,8,7,4,6,4,8,8,11,10,6,12,5,4,5,13,6,12},
		{13,0,7,1,5,15,13,2,8,3,2,3,12,12,14,12,2,1,12,14,15},
		{10,7,0,4,8,12,7,11,14,8,11,6,5,11,9,3,14,5,15,4,14},
		{14,1,4,0,2,3,11,2,15,5,15,10,13,4,13,11,8,2,7,3,4},
		{8,5,8,2,0,4,5,3,3,9,12,13,8,15,2,15,14,3,12,8,15},
		{7,15,12,3,4,0,15,5,15,11,15,13,10,13,3,12,7,1,5,8,3},
		{4,13,7,11,5,15,0,15,7,7,13,7,7,6,8,13,3,5,6,10,9},
		{6,2,11,2,3,5,15,0,8,13,5,7,10,1,3,1,9,14,13,3,7},
		{4,8,14,15,3,15,7,8,0,14,11,3,2,8,9,5,6,15,14,14,13},
		{8,3,8,5,9,11,7,13,14,0,12,3,12,7,2,14,1,8,9,2,6},
		{8,2,11,15,12,15,13,5,11,12,0,3,2,10,6,13,13,11,5,3,9},
		{11,3,6,10,13,13,7,7,3,3,3,0,11,13,8,9,2,1,6,10,9},
		{10,12,5,13,8,10,7,10,2,12,2,11,0,13,13,2,12,13,14,11,15},
		{6,12,11,4,15,13,6,1,8,7,10,13,13,0,15,15,12,11,9,4,14},
		{12,14,9,13,2,3,8,3,9,2,6,8,13,15,0,8,11,9,11,9,12},
		{5,12,3,11,15,12,13,1,5,14,13,9,2,15,8,0,11,5,11,14,13},
		{4,2,14,8,14,7,3,9,6,1,13,2,12,12,11,11,0,2,12,10,7},
		{5,1,5,2,3,1,5,14,15,8,11,1,13,11,9,5,2,0,9,11,5},
		{13,12,15,7,12,5,6,13,14,9,5,6,14,9,11,11,12,9,0,6,9},
		{6,14,4,3,8,8,10,3,14,2,3,10,11,4,9,14,10,11,6,0,2},
		{12,15,14,4,15,3,9,7,13,6,9,9,15,14,12,13,7,5,9,2,0}
	};
	for (int i = 0; i < C; i++)
	{
		for (int j = 0; j < C; j++)
		{
			Dist[i][j] = d2[i][j];
		}
	}
}

/*《函数ValueInit：初始化各变量数组》*/
void ValueInit()
{

	for (int i = 0; i < C; i++)			//初始化 Eta[n][n]，地点间距离的倒数
	{
		for (int j = 0; j < C; j++)
		{
			Eta[i][j] = 1.0 / Dist[i][j];
		}
	}
	for (int i = 0; i < C; i++)			//初始化 Tau[n][n]，两两节点间的初始信息素浓度
	{
		for (int j = 0; j < C; j++)
		{
			Tau[i][j] = 1.0;
		}
	}
	for (int i = 0; i < C; i++)			//初始化 DeltaEta[n][n],两两节点间信息素的变化量
	{
		for (int j = 0; j < C; j++)
		{
			DeltaTau[i][j] = 0;
		}
	}
	for (int i = 0; i < Ants; i++)		//初始化 Tabu[m][n]，禁忌表
	{
		for (int j = 0; j < C; j++)
		{
			Tabu[i][j] = 0;
		}
	}
}
/*《函数Random_Num：生成lower和uper之间的一个double类型随机数》*/
double Random_Num(double lower, double uper)
{
	return  (rand() / (double)RAND_MAX) * (uper - lower) + lower;
}

int all_citys[C];//生成一个待访问地点序列
int main()
{


	bool x_h = true;
	string ljd;
	string qi_dian;
	bool cw;
	while (x_h)
	{
		for (int i = 0; i < 2; i++)
		{
			cout << i << "," << qi_dian_b[i] << "  ";
		}
		cout << "请输入起点：";
		cin >> ljd;
		cw = true;
		int cz = find(qi_dian_b, qi_dian_b + 2, ljd) - qi_dian_b;
		if (cz != 2)
		{
			di_dian[0] = ljd;
			int wz = find(di_dian, di_dian + 22, ljd) - di_dian;
			lu_jing_dian[lu_jing_dian_shu - 1] = wz;
			qi_dian = ljd;
			lu_jing_dian_shu += 1;
			if (wz == 1)
			{
				hb();
			}
			break;
			cw = false;
		}

		if (cw)
		{
			lu_jing_dian_shu -= 1;
			cout << "起点未知或不可用，请重新输入" << endl;
		}

	}


	while (x_h)
	{
		lu_jing_dian_shu += 1;

		for (int i = 0; i < C; i++)
		{
			cout << i << "," << di_dian[i] << "  ";
		}
		cout << endl << "请输入路径点" << lu_jing_dian_shu << "：";
		cin >> ljd;

		if (ljd == qi_dian)
		{
			cout << "所选路径点为：";
			for (int i = 0; i < lu_jing_dian_shu - 1; i++)
			{
				cout << di_dian[lu_jing_dian[i]] << "  ";
			}
			cout << endl;
			lu_jing_dian_shu -= 1;
			break;
		}

		cw = true;
		int wz = find(di_dian, di_dian + 22, ljd) - di_dian;
		cout << wz << endl;
		if (wz != C)
		{
			int cf = find(lu_jing_dian, lu_jing_dian + C, wz) - lu_jing_dian;
			cout << cf;
			if (cf == C)
			{
				lu_jing_dian[lu_jing_dian_shu - 1] = wz;
				cw = false;
			}
		}

		if (cw)
		{
			lu_jing_dian_shu -= 1;
			cout << "路径点未知或重复，请重新输入" << endl;
		}
	}



	int i, j;
	/*Step1：调用函数《ValueInit》:初始化各变量数组*/
	ValueInit();

	int it = 0;//迭代参数
	while (it < I)
	{
		//for (i = 0; i < c; i++)
		//{
		//	all_citys[i] = i + 1;				//生成一个待访问地点全序列	
		//}
		for (i = 0; i < C; i++)
		{
			all_citys[i] = lu_jing_dian[i];
		}
		

		/*Step2：将 Ants 只蚂蚁随机平均放到 C 个地点上*/
		//①一维向量，连续存放几轮对各地点安放蚂蚁的结果
		vector<int> temp;
		//②ceil()函数向上取整；floor()函数向下取整；round()四舍五入
		for (int i = 0; i < ceil((double)Ants / (double)C); i++)
		{
			for (int j = 0; j < C; j++)
				temp.push_back(j);
		}
		//③打乱temp数组中元素的次序
		//random_shuffle(temp.begin(), temp.end());
		
		//④更新每只蚂蚁的禁忌表，实际只用到temp前m个节点的信息
		//for (int i = 0; i < Ants; i++)
		//{
		//	Tabu[i][0] = temp[i];
		//}

		
		for (int i = 0; i < Ants; i++)
		{
			int rand_city = 0;			//随机为当前蚂蚁选择一个出发地点索引
			Tabu[i][0] = all_citys[rand_city];	//同时更新每只蚂蚁的禁忌表
		}
		
		/*Step3：Ants只蚂蚁按概率函数选择下一座地点进行访问，并完成各自的周游*/
		//①禁忌表已经有第一个出发地点，所以从第二个点开始概率选择访问（基于信息素和启发因子）
		for (int j = 1; j < C; j++)
		{
			//②Ants只蚂蚁并行
			for (int i = 0; i < Ants; i++)
			{
				vector<int> visited;		//第i只蚂蚁，已访问过的地点，集合
				vector<int> J;				//第i只蚂蚁，待访问的地点，集合
				vector<double> V;				//第i只蚂蚁，待访问地点的竞争值，的集合
				vector<double> P;			//第i只蚂蚁，待访问地点的被选概率，的集合

				double Vsum = 0.0;			//竞争值和
				double Psum = 0.0;			//概率和
				double rand_rate = 0.0;		//随机概率
				double choose = 0.0;		//轮盘赌累加概率（用于和rand_rate比较）
				int to_visit = 0;			//下一个要去的地点

				//③visited，结算当前蚂蚁已访问节点
				for (int k = 0; k < j; k++)			//j控制当前已访问节点个数（包含出发点）
				{
					visited.push_back(Tabu[i][k]);	//把当前蚂蚁的禁忌表拿来赋值
				}

				//④J，添加待访问节点
				for (int k = 0; k < C; k++)
				{
					//若在visited中没有找到某节点all_citys[k]，则说明还未被访问
					if (find(visited.begin(), visited.end(), all_citys[k]) == visited.end())
					{
						J.push_back(all_citys[k]);			//J初始化，添加到待访问地点集合
						P.push_back(0.0);					//P初始化，待访问地点被选概率占坑
						V.push_back(0.0);					//V初始化，待访问地点竞争值占坑
					}
				}

				//⑤计算各待访地点的竞争值（已访问节点中最后一个地点—>各待访问各节点）
				for (int k = 0; k < V.size(); k++)
				{
					V[k] = pow(Tau[visited.back()][J[k]], Alpha) * pow(Eta[visited.back()][J[k]], Beta);
					Vsum += V[k];
				}

				//⑥计算去下一座地点的各种概率（已访问节点中最后一个地点—>各待访问各节点）
				for (int k = 0; k < P.size(); k++)
				{
					P[k] = V[k] / Vsum;
					Psum += P[k];
				}

				//⑦使用轮盘赌算法，挑选下一座要去的地点
				rand_rate = Random_Num(0.0, Psum);		//随机概率		
				for (int k = 0; k < P.size(); k++)
				{
					choose += P[k];						//概率小累计，用于与随机概率比较
					if (choose > rand_rate)
					{
						to_visit = J[k];		//根据概率选择下一访问地点
						break;
					}
				}

				//⑧更新禁忌表（把刚选出的下一访问节点及时加入）
				Tabu[i][j] = to_visit;

			}
		}

		//Step4：计算本次迭代各蚂蚁的旅行数据（距离、时间等）
		//①记录本代每只蚂蚁走的总路程
		double L[Ants];
		for (int i = 0; i < Ants; i++)
		{
			L[i] = 0.0;						//初始化
		}
		for (int i = 0; i < Ants; i++)
		{
			for (int j = 0; j < C - 1; j++)
			{
				//禁忌表中就是最终的访问路径（但禁忌表从1开始，而Dist[][]从0开始存储，所以要减一）
				L[i] += Dist[Tabu[i][j] - 1][Tabu[i][j + 1] - 1];
			}
			L[i] += Dist[Tabu[i][0] - 1][Tabu[i][C - 1] - 1];
		}
		//②计算当前迭代所有蚂蚁路径中的最小蚂蚁路径
		double min_value = L[0];					//声明求本代所有蚂蚁行走距离最小值的临时变量
		double sum_value = L[0];					//声明求本代所有蚂蚁行走距离总值的临时变量
		int min_index = 0;							//记录本代所有蚂蚁行走距离最小值的下标
		for (int i = 1; i < Ants; i++)
		{
			sum_value += L[i];
			if (L[i] < min_value)
			{
				min_value = L[i];
				min_index = i;
			}
		}
		//③将本次迭代的最小路径值、平均路径值、最短路径数据存在全局迭代数组中
		L_best[it] = min_value;						//每代中路径的最短长度
		L_ave[it] = sum_value / Ants;				//每代中路径的平均长度
		for (int i = 0; i < C; i++)
		{
			R_best[it][i] = Tabu[min_index][i];		//记录每代最短的路径数据
		}
		//④打印各代距离信息
		//cout << it << ": L_best is " << L_best[it] << "    " << "L_ave is " << L_ave[it] << endl;

		//⑤迭代继续
		it++;

		//Step5：更新两两节点间的信息素，为下一次迭代的概率选择下一访问节点作准备
		//①全部更新：累加所有蚂蚁，所有路径中，涉及到相同的两两节点间的信息素增量
		for (int i = 0; i < Ants; i++)
		{
			for (int j = 0; j < C - 1; j++)
			{
				//（禁忌表从1开始，而DeltaTau[][]从0开始存储，所以要减一）
				DeltaTau[Tabu[i][j] - 1][Tabu[i][j + 1] - 1] += Q / L[i];
			}
			//蚁周算法，使用全局距离作更新
			DeltaTau[Tabu[i][C - 1] - 1][Tabu[i][0] - 1] += Q / L[i];
		}
		//②基于以上两两节点间的信息素增量DeltaTau[][]，对Tau[][]作更新
		for (int i = 0; i < C; i++)
		{
			for (int j = 0; j < C; j++)
			{
				//考虑信息素挥发+信息素增量，更新信息素
				Tau[i][j] = (1 - Rho) * Tau[i][j] + DeltaTau[i][j];
			}
		}
		//③禁忌表清零,进入下一次迭代!
		for (int i = 0; i < Ants; i++)
		{
			for (int j = 0; j < C; j++)
			{
				Tabu[i][j] = 0;
			}
		}

	}

	//Step6：展示最终所有迭代中最好的结果
	double min_L = L_best[0];			//所有迭代中最短距离
	int min_L_index = 0;				//所有迭代中最优路径的迭代下标
	int Shortest_Route[C];				//所有迭代中的最优路径
	for (int i = 0; i < it; i++)
	{
		if (L_best[i] < min_L)
		{
			min_L = L_best[i];
			min_L_index = i;
		}
	}
	cout << endl << endl;
	//cout << "The length of the shortest route is :  " << min_L << endl;
	//cout << "The number of iteration is :  " << min_L_index << endl;
	//cout << "The Shortest route is： " << endl;// << "[start]";

	//for (int i = 0; i < C; i++)		//所有迭代中的最优路径
	//{
	//	Shortest_Route[i] = R_best[min_L_index][i];
	//	cout << " -> " << Shortest_Route[i];
	//}
	for (int i = 0; i < lu_jing_dian_shu; i++)
	{
		Shortest_Route[i] = R_best[min_L_index][i];
	}
	int xhks = find(Shortest_Route, Shortest_Route + C, 0) - Shortest_Route;
	cout << di_dian[Shortest_Route[xhks]];
	for (int i = 1; i < lu_jing_dian_shu; i++)		//所有迭代中的最优路径
	{
		if (i + xhks < lu_jing_dian_shu)
		{
			cout << " -> " << di_dian[Shortest_Route[i + xhks]];
		}
		else
		{
			cout << " -> " << di_dian[Shortest_Route[i + xhks - lu_jing_dian_shu]];
		}

	}
	cout << " -> " << di_dian[Shortest_Route[xhks]] << endl;

	system("pause");
	return 0;
}
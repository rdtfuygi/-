#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
using namespace std;

#define Ants 100			//���ϵĸ���
#define C 21				//�ص������
#define I 50				//����������
#define Alpha 1				//������Ϣ����Ҫ�̶ȵĲ���
#define Beta 5				//��������ʽ������Ҫ�̶ȵĲ���
#define Rho 0.1				//��Ϣ������ϵ��
#define Q 100				//��Ϣ������ǿ��ϵ��

string qi_dian_b[2] = { "��ʳ","��ʳ" };
string di_dian[C] = { "��ʳ","��1","��2","��3","��4","��5","��6","��7","��8","��9","��10" "��1","��2","��3","��4","��5","��6","��7","��8","��9","��10","��11" };
int lu_jing_dian_shu = 0;
int lu_jing_dian[C];


int Tabu[Ants][C];			//���ɱ��洢�߹���·��
int Dist[C][C] =		//��ʾ�����ص�����,ֱ���ò���ʱ���ʾ
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
double Eta[C][C];			//��ʾ����ʽ���ӣ�ΪDist[][]�о���ĵ���
double Tau[C][C];			//��ʾ�����ڵ����Ϣ��Ũ��
double DeltaTau[C][C];		//��ʾ�����ڵ����Ϣ�صı仯��
double L_best[I];			//�洢ÿ�ε�����·������̳���
double L_ave[I];			//�洢ÿ�ε�����·����ƽ������
int R_best[I][C];			//�洢ÿ�ε��������·��

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

/*������ValueInit����ʼ�����������顷*/
void ValueInit()
{

	for (int i = 0; i < C; i++)			//��ʼ�� Eta[n][n]���ص�����ĵ���
	{
		for (int j = 0; j < C; j++)
		{
			Eta[i][j] = 1.0 / Dist[i][j];
		}
	}
	for (int i = 0; i < C; i++)			//��ʼ�� Tau[n][n]�������ڵ��ĳ�ʼ��Ϣ��Ũ��
	{
		for (int j = 0; j < C; j++)
		{
			Tau[i][j] = 1.0;
		}
	}
	for (int i = 0; i < C; i++)			//��ʼ�� DeltaEta[n][n],�����ڵ����Ϣ�صı仯��
	{
		for (int j = 0; j < C; j++)
		{
			DeltaTau[i][j] = 0;
		}
	}
	for (int i = 0; i < Ants; i++)		//��ʼ�� Tabu[m][n]�����ɱ�
	{
		for (int j = 0; j < C; j++)
		{
			Tabu[i][j] = 0;
		}
	}
}
/*������Random_Num������lower��uper֮���һ��double�����������*/
double Random_Num(double lower, double uper)
{
	return  (rand() / (double)RAND_MAX) * (uper - lower) + lower;
}

int all_citys[C];//����һ�������ʵص�����
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
		cout << "��������㣺";
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
			cout << "���δ֪�򲻿��ã�����������" << endl;
		}

	}


	while (x_h)
	{
		lu_jing_dian_shu += 1;

		for (int i = 0; i < C; i++)
		{
			cout << i << "," << di_dian[i] << "  ";
		}
		cout << endl << "������·����" << lu_jing_dian_shu << "��";
		cin >> ljd;

		if (ljd == qi_dian)
		{
			cout << "��ѡ·����Ϊ��";
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
			cout << "·����δ֪���ظ�������������" << endl;
		}
	}



	int i, j;
	/*Step1�����ú�����ValueInit��:��ʼ������������*/
	ValueInit();

	int it = 0;//��������
	while (it < I)
	{
		//for (i = 0; i < c; i++)
		//{
		//	all_citys[i] = i + 1;				//����һ�������ʵص�ȫ����	
		//}
		for (i = 0; i < C; i++)
		{
			all_citys[i] = lu_jing_dian[i];
		}
		

		/*Step2���� Ants ֻ�������ƽ���ŵ� C ���ص���*/
		//��һά������������ż��ֶԸ��ص㰲�����ϵĽ��
		vector<int> temp;
		//��ceil()��������ȡ����floor()��������ȡ����round()��������
		for (int i = 0; i < ceil((double)Ants / (double)C); i++)
		{
			for (int j = 0; j < C; j++)
				temp.push_back(j);
		}
		//�۴���temp������Ԫ�صĴ���
		//random_shuffle(temp.begin(), temp.end());
		
		//�ܸ���ÿֻ���ϵĽ��ɱ�ʵ��ֻ�õ�tempǰm���ڵ����Ϣ
		//for (int i = 0; i < Ants; i++)
		//{
		//	Tabu[i][0] = temp[i];
		//}

		
		for (int i = 0; i < Ants; i++)
		{
			int rand_city = 0;			//���Ϊ��ǰ����ѡ��һ�������ص�����
			Tabu[i][0] = all_citys[rand_city];	//ͬʱ����ÿֻ���ϵĽ��ɱ�
		}
		
		/*Step3��Antsֻ���ϰ����ʺ���ѡ����һ���ص���з��ʣ�����ɸ��Ե�����*/
		//�ٽ��ɱ��Ѿ��е�һ�������ص㣬���Դӵڶ����㿪ʼ����ѡ����ʣ�������Ϣ�غ��������ӣ�
		for (int j = 1; j < C; j++)
		{
			//��Antsֻ���ϲ���
			for (int i = 0; i < Ants; i++)
			{
				vector<int> visited;		//��iֻ���ϣ��ѷ��ʹ��ĵص㣬����
				vector<int> J;				//��iֻ���ϣ������ʵĵص㣬����
				vector<double> V;				//��iֻ���ϣ������ʵص�ľ���ֵ���ļ���
				vector<double> P;			//��iֻ���ϣ������ʵص�ı�ѡ���ʣ��ļ���

				double Vsum = 0.0;			//����ֵ��
				double Psum = 0.0;			//���ʺ�
				double rand_rate = 0.0;		//�������
				double choose = 0.0;		//���̶��ۼӸ��ʣ����ں�rand_rate�Ƚϣ�
				int to_visit = 0;			//��һ��Ҫȥ�ĵص�

				//��visited�����㵱ǰ�����ѷ��ʽڵ�
				for (int k = 0; k < j; k++)			//j���Ƶ�ǰ�ѷ��ʽڵ���������������㣩
				{
					visited.push_back(Tabu[i][k]);	//�ѵ�ǰ���ϵĽ��ɱ�������ֵ
				}

				//��J����Ӵ����ʽڵ�
				for (int k = 0; k < C; k++)
				{
					//����visited��û���ҵ�ĳ�ڵ�all_citys[k]����˵����δ������
					if (find(visited.begin(), visited.end(), all_citys[k]) == visited.end())
					{
						J.push_back(all_citys[k]);			//J��ʼ������ӵ������ʵص㼯��
						P.push_back(0.0);					//P��ʼ���������ʵص㱻ѡ����ռ��
						V.push_back(0.0);					//V��ʼ���������ʵص㾺��ֵռ��
					}
				}

				//�ݼ�������õص�ľ���ֵ���ѷ��ʽڵ������һ���ص㡪>�������ʸ��ڵ㣩
				for (int k = 0; k < V.size(); k++)
				{
					V[k] = pow(Tau[visited.back()][J[k]], Alpha) * pow(Eta[visited.back()][J[k]], Beta);
					Vsum += V[k];
				}

				//�޼���ȥ��һ���ص�ĸ��ָ��ʣ��ѷ��ʽڵ������һ���ص㡪>�������ʸ��ڵ㣩
				for (int k = 0; k < P.size(); k++)
				{
					P[k] = V[k] / Vsum;
					Psum += P[k];
				}

				//��ʹ�����̶��㷨����ѡ��һ��Ҫȥ�ĵص�
				rand_rate = Random_Num(0.0, Psum);		//�������		
				for (int k = 0; k < P.size(); k++)
				{
					choose += P[k];						//����С�ۼƣ�������������ʱȽ�
					if (choose > rand_rate)
					{
						to_visit = J[k];		//���ݸ���ѡ����һ���ʵص�
						break;
					}
				}

				//����½��ɱ��Ѹ�ѡ������һ���ʽڵ㼰ʱ���룩
				Tabu[i][j] = to_visit;

			}
		}

		//Step4�����㱾�ε��������ϵ��������ݣ����롢ʱ��ȣ�
		//�ټ�¼����ÿֻ�����ߵ���·��
		double L[Ants];
		for (int i = 0; i < Ants; i++)
		{
			L[i] = 0.0;						//��ʼ��
		}
		for (int i = 0; i < Ants; i++)
		{
			for (int j = 0; j < C - 1; j++)
			{
				//���ɱ��о������յķ���·���������ɱ��1��ʼ����Dist[][]��0��ʼ�洢������Ҫ��һ��
				L[i] += Dist[Tabu[i][j] - 1][Tabu[i][j + 1] - 1];
			}
			L[i] += Dist[Tabu[i][0] - 1][Tabu[i][C - 1] - 1];
		}
		//�ڼ��㵱ǰ������������·���е���С����·��
		double min_value = L[0];					//�����󱾴������������߾�����Сֵ����ʱ����
		double sum_value = L[0];					//�����󱾴������������߾�����ֵ����ʱ����
		int min_index = 0;							//��¼���������������߾�����Сֵ���±�
		for (int i = 1; i < Ants; i++)
		{
			sum_value += L[i];
			if (L[i] < min_value)
			{
				min_value = L[i];
				min_index = i;
			}
		}
		//�۽����ε�������С·��ֵ��ƽ��·��ֵ�����·�����ݴ���ȫ�ֵ���������
		L_best[it] = min_value;						//ÿ����·������̳���
		L_ave[it] = sum_value / Ants;				//ÿ����·����ƽ������
		for (int i = 0; i < C; i++)
		{
			R_best[it][i] = Tabu[min_index][i];		//��¼ÿ����̵�·������
		}
		//�ܴ�ӡ����������Ϣ
		//cout << it << ": L_best is " << L_best[it] << "    " << "L_ave is " << L_ave[it] << endl;

		//�ݵ�������
		it++;

		//Step5�����������ڵ�����Ϣ�أ�Ϊ��һ�ε����ĸ���ѡ����һ���ʽڵ���׼��
		//��ȫ�����£��ۼ��������ϣ�����·���У��漰����ͬ�������ڵ�����Ϣ������
		for (int i = 0; i < Ants; i++)
		{
			for (int j = 0; j < C - 1; j++)
			{
				//�����ɱ��1��ʼ����DeltaTau[][]��0��ʼ�洢������Ҫ��һ��
				DeltaTau[Tabu[i][j] - 1][Tabu[i][j + 1] - 1] += Q / L[i];
			}
			//�����㷨��ʹ��ȫ�־���������
			DeltaTau[Tabu[i][C - 1] - 1][Tabu[i][0] - 1] += Q / L[i];
		}
		//�ڻ������������ڵ�����Ϣ������DeltaTau[][]����Tau[][]������
		for (int i = 0; i < C; i++)
		{
			for (int j = 0; j < C; j++)
			{
				//������Ϣ�ػӷ�+��Ϣ��������������Ϣ��
				Tau[i][j] = (1 - Rho) * Tau[i][j] + DeltaTau[i][j];
			}
		}
		//�۽��ɱ�����,������һ�ε���!
		for (int i = 0; i < Ants; i++)
		{
			for (int j = 0; j < C; j++)
			{
				Tabu[i][j] = 0;
			}
		}

	}

	//Step6��չʾ�������е�������õĽ��
	double min_L = L_best[0];			//���е�������̾���
	int min_L_index = 0;				//���е���������·���ĵ����±�
	int Shortest_Route[C];				//���е����е�����·��
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
	//cout << "The Shortest route is�� " << endl;// << "[start]";

	//for (int i = 0; i < C; i++)		//���е����е�����·��
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
	for (int i = 1; i < lu_jing_dian_shu; i++)		//���е����е�����·��
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
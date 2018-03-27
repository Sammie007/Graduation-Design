#ifndef _D
#define _D

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctime>
#include <sys/timeb.h>
#include <cerrno>
#include <unistd.h>
#include <csignal>
#include <cmath>
#include <algorithm>

using namespace std;

const int     MAX_LINE_LEN           = 50;
const int     MAX_SENSOR_NUM         = 10000;
const int     DEF_CHROMOSOME_NUM     = 20;             //染色体规模,默认值为20
const int     DEF_MAX_GENERATIONS    = 100;          //遗传算法迭代的代数
const double  DEF_PROBABILITY_MUTATE = 0.1;       //默认变异概率
const double  DEF_PROBABILITY_CROSSOVER = 0.85;    //默认交叉概率
const double  DEF_MOBILE_PERCENT = 0.30;           //移动节点比例
const int     DEF_GRID_NUM = 100;
const int     MAX_MOBILE_NUM = 1000000;
//const int     DEF_SENSOR_NUM = 60;              //缺省的节点个数(包含移动节点和固定节点)

int     chromosome_num  = DEF_CHROMOSOME_NUM;
int     max_generations = DEF_MAX_GENERATIONS;
double  probability_mutate = DEF_PROBABILITY_MUTATE;
double  probability_crossover = DEF_PROBABILITY_CROSSOVER;
double  mobile_percent = DEF_MOBILE_PERCENT;
//int     point_radius = 3;
int     grid = 28;
int     field[4];    //0和1表示x轴的范围(默认190~470)，2和3表示y轴的范围(默认84~364)
int     row_grid_num, col_grid_num;
int     grid_num;     //平面上单元格的数目
bool    cover_flag[DEF_GRID_NUM];       //两个状态，1-覆盖，0-未覆盖
bool    static_cover[DEF_GRID_NUM];
int     placement_unit;
int     max_move_dist;
//double  sigma = 0.1;
double  borderlength;
int     sensor_num;
double  transmission_range;
double  sense_range;
double  ideal_length;
double  max_coverage;
double  move_d;


bool    is_mobile[MAX_SENSOR_NUM];       //是否是移动节点
int     mobile_num;
int     mobile_pos[MAX_MOBILE_NUM];       //记录移动节点在node中的位置

const double PI = 3.1415926;
const double E = 2.7182818;

struct Node{
   int x, y;   //x,y坐标
};

struct Chromosomes{
   vector<int> x;
   vector<int> y;
   int label;   //属性序号
};

struct chromo_attribute{    //染色体的属性
   double coverage;    //覆盖率
   double init_c;      //覆盖率副本
   int label;          //染色体序号
   bool is_select = 0;      //染色体是否被选中
   double move_d = 0.0;       //移动距离
};

Node  node[MAX_SENSOR_NUM];                    //传感器节点
Chromosomes chromosome[DEF_CHROMOSOME_NUM];        //染色体
chromo_attribute des[DEF_CHROMOSOME_NUM];

int  read_file(const char * const filename);
void produce_chrom();
double compute_coverage();
void cover(Node node2);
void static_covere();
void evaluate();
void select();
void select2();
void cross();
void variation(int generation);
double evaluate_single(Chromosomes chro);
void single_point_crossover();
void cover_first();
#endif // _D


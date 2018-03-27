#include "mainwindow.h"
#include "1.h"
#include <QApplication>
#include <QMessageBox>
#include <QFileDialog>
#include <QDebug>
#include <QTime>
#include <QDir>
#include <QFile>
#include <iostream>
#include <vector>

bool cmp1(const chromo_attribute &a, const chromo_attribute &b)
{
    return a.coverage > b.coverage;
}

bool cmp2(const chromo_attribute &a, const chromo_attribute &b)
{
    if(a.is_select != b.is_select) return a.is_select > b.is_select;
    else return a.coverage > b.coverage ;
}

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}

int read_file(const char * const filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Fail to open file %s, %s.\n", filename, strerror(errno));
        return 0;
    }
    printf("Open file %s OK.\n", filename);

    unsigned int cnt = 0, line_num = 0;
    char line[MAX_LINE_LEN + 2];
    char temp[10]; for(int i=0;i<10;i++) temp[i] = 0;
    int node_label;

    line[0] = 0;
    if (fgets(line, MAX_LINE_LEN + 2, fp) == NULL || line[0] == 0)  return -1;
    else sscanf(line, "%d", &line_num);

    while ((cnt < line_num) && !feof(fp))
    {
        line[0] = 0;
        if (fgets(line, MAX_LINE_LEN + 2, fp) == NULL)  continue;
        if (line[0] == 0)  continue;

        int i,j;
        for(i=0;i<strlen(line) && line[i] != ' ';i++)
             temp[i] = line[i];
        sscanf(temp, "%d", &node_label);    //找到结点序号
        for(j=0;j<10;j++) temp[j] = 0;

        while(line[i] == ' ') i++;
        for(j=0;i<strlen(line) && line[i] != ' ';i++,j++)
             temp[j] = line[i];
        sscanf(temp, "%d", &node[node_label-1].x);    //x坐标
        for(j=0;j<10;j++) temp[j] = 0;

        while(line[i] == ' ') i++;
        for(j=0;i<strlen(line) && line[i] != ' ';i++,j++)
             temp[j] = line[i];
        sscanf(temp, "%d", &node[node_label-1].y);     //y坐标
        for(j=0;j<10;j++) temp[j] = 0;

        cnt++;
    }
    fclose(fp);
    printf("There are %d lines in file %s.\n", cnt, filename);

    return line_num;
}

void variation(int generation)
{
    //有chromosome_num条染色体，每条有mobile_num个移动节点，循环一次变异一个移动节点
    //故有chromosome_num * mobile_num次操作
     for(int chrom=0;chrom<chromosome_num;chrom++){
        double p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
        if(p > DEF_PROBABILITY_MUTATE)    continue;  //变异概率

        for(int mobile = 0;mobile < mobile_num;mobile++){
             p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
             if(p > DEF_PROBABILITY_MUTATE)    continue;  //变异概率
//cout<<"variation"<<endl;
             //int  chrom   = i;     //0到chromosome_num - 1, 即选中一条染色体
             //int mobile   = i % mobile_num;    //0到mobile_num - 1, 即选中染色体中的一个传感器
//cout<<chrom<<' '<<mobile<<endl;
             double angle = (rand()%(628 - 0 + 1) + 0) / 100.0;
             //随机产生的角度, 因一圈为2π，即2 * 3.14 = 6.28
             double r     = (1 - generation / max_generations) * ideal_length;
             int    x     = chromosome[chrom].x[mobile] + (int)(r * cos(angle) + 0.5);  //+0.5再取整是为了四舍五入
             int    y     = chromosome[chrom].y[mobile] - (int)(r * sin(angle) + 0.5);
             x = x < field[0] ? field[0] : x;
             x = x > field[1] ? field[1] : x;
             y = y < field[2] ? field[2] : y;
             y = y > field[3] ? field[3] : y;
//cout<<chromosome[chrom].x[mobile]<<' '<<x<<"    "<<chromosome[chrom].y[mobile]<<' '<<y<<endl;
             chromosome[chrom].x[mobile] = x;
             chromosome[chrom].y[mobile] = y;

        }
     }
}


void cross()
{
    for(int i=0;i<chromosome_num / 2;i++){   //循环chromosome_num / 2次
        double p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
        if(p > DEF_PROBABILITY_CROSSOVER) continue;    //交叉概率
        Chromosomes temp1, temp2;
        int one = des[2 * i].label;        //找出排序后的染色体的真实位置
        int two = des[2 * i + 1].label;
//cout<<one<<' '<<two<<' ';
        int m=0;
        //m = rand() % (mobile_num + 1);
        for(;m<mobile_num;m++){
            double p1 = (double)(rand() * 1.0 /(RAND_MAX+1.0));      //产生从0到1（不包括1）的随机数
            double p2 = (double)(rand() * 1.0 /(RAND_MAX+1.0));
            int t1 = chromosome[one].x[m];     //暂存
            int t2 = chromosome[two].x[m];     //暂存
            int t3 = chromosome[one].y[m];     //暂存
            int t4 = chromosome[two].y[m];     //暂存
            p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
            if(p < DEF_PROBABILITY_CROSSOVER)
                temp1.x.push_back(round( p1   * t1 + (1-p1)* t2));
            else temp1.x.push_back(chromosome[one].x[m]);

            p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
            if(p < DEF_PROBABILITY_CROSSOVER)
                temp2.x.push_back(round((1-p1)* t1 +  p1   * t2));
            else temp2.x.push_back(chromosome[two].x[m]);

            p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
            if(p < DEF_PROBABILITY_CROSSOVER)
                temp1.y.push_back(round( p2   * t3 + (1-p2)* t4));
            else temp1.y.push_back(chromosome[one].y[m]);

            p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
            if(p < DEF_PROBABILITY_CROSSOVER)
                temp2.y.push_back(round((1-p2)* t3 +  p2   * t4));
            else temp2.y.push_back(chromosome[two].y[m]);

//cout<<"before: "<<chromosome[one].x[m]<<' '<<chromosome[two].x[m]<<' '<<chromosome[one].y[m]<<' '<<chromosome[two].y[m]<<endl;
            //chromosome[one].x[m] = t_x1;
            //chromosome[two].x[m] = t_x2;
            //chromosome[one].y[m] = t_y1;
            //chromosome[two].y[m] = t_y2;
//cout<<"after: "<<chromosome[one].x[m]<<' '<<chromosome[two].x[m]<<' '<<chromosome[one].y[m]<<' '<<chromosome[two].y[m]<<endl;cout<<endl;
        }
        double t1 = evaluate_single(temp1);
        double t2 = evaluate_single(temp2);
        if(t1 > des[2 * i].coverage){
//cout<<"Yes1 "<<t1<<' '<<des[2 * i].coverage<<endl;
            for(int m=0;m<mobile_num;m++){
                chromosome[one].x[m] = temp1.x[m];
                chromosome[one].y[m] = temp1.y[m];
            }
        }
        //else cout<<"No1"<<endl;
        if( t2 > des[2 * i + 1].coverage){
//cout<<"Yes2 "<<t2<<' '<<des[2 * i + 1].coverage<<endl;
            for(int m=0;m<mobile_num;m++){
                chromosome[two].x[m] = temp2.x[m];
                chromosome[two].y[m] = temp2.y[m];
            }
        }
        //else cout<<"No2"<<endl;
    }

}

void single_point_crossover()
{
    //one和two为参与交叉的两条染色体
    int one = -1, two = -1;
    for(int i=0;i<chromosome_num;i++){
        double p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
        if(p > DEF_PROBABILITY_CROSSOVER) continue;
        Chromosomes temp1, temp2;
        if(one == -1) one = i;
        else if(two == -1){
            two = i;
//cout<<one<<' '<<two<<endl;
            //int pos = rand() % (mobile_num - 1);    //pos为0到mobile_num - 1
            for(int m=0;m<mobile_num;m++){
//cout<<chromosome[one].x[m]<<' '<<chromosome[two].x[m]<<' ';
                 p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
                 if(p < DEF_PROBABILITY_CROSSOVER)
                        temp1.x.push_back(chromosome[two].x[m]);
                 else temp1.x.push_back(chromosome[one].x[m]);

                 p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
                 if(p < DEF_PROBABILITY_CROSSOVER)
                        temp1.y.push_back(chromosome[two].y[m]);
                 else temp1.y.push_back(chromosome[one].y[m]);

                 p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
                 if(p < DEF_PROBABILITY_CROSSOVER)
                        temp2.x.push_back(chromosome[one].x[m]);
                 else temp2.x.push_back(chromosome[two].x[m]);

                 p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
                 if(p < DEF_PROBABILITY_CROSSOVER)
                        temp2.y.push_back(chromosome[one].y[m]);
                 else temp2.y.push_back(chromosome[two].y[m]);
                 //swap(chromosome[two].x[m] , chromosome[one].x[m]);
                 //swap(chromosome[two].y[m] , chromosome[one].y[m]);
//cout<<chromosome[one].x[m]<<' '<<chromosome[two].x[m]<<endl;
            }
        double t1 = evaluate_single(temp1);
        double t2 = evaluate_single(temp2);
        //此步非常重要，交叉后比原来覆盖率高才执行，否则保留优质个体
        if(t1 > des[chromosome[one].label].coverage){
//cout<<"Yes1 "<<t1<<' '<<des[2 * i].coverage<<endl;
            for(int m=0;m<mobile_num;m++){
                chromosome[one].x[m] = temp1.x[m];
                chromosome[one].y[m] = temp1.y[m];
            }
        }
        //else cout<<"No1"<<endl;
        if( t2 > des[chromosome[two].label].coverage){
//cout<<"Yes2 "<<t2<<' '<<des[2 * i + 1].coverage<<endl;
            for(int m=0;m<mobile_num;m++){
                chromosome[two].x[m] = temp2.x[m];
                chromosome[two].y[m] = temp2.y[m];
            }
        }
           one = two = -1;
        }
    }
}

void select()
{
     double sum = 0, select_p, accu = 0;   //accu为累积选中概率
     for(int i=0;i<chromosome_num;i++) sum += des[i].coverage;
     des[0].is_select = true;      //最大适应值的总被选中
     for(int i=1;i<chromosome_num;i++){
        select_p = des[i].coverage / sum;
        accu += select_p;
        double p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
        //选择概率，等于适应值大小除以适应值总和
        if( p > accu){
//cout<<p<<' '<<accu<<endl;
            des[i].is_select = true;
        }
     }
     sort(des, des + chromosome_num, cmp2);
}

void select2()
{
/*
for(int i=0;i<chromosome_num;i++){
    for(int j=0;j<chromosome[i].x.size();j++){
        cout<<chromosome[i].x[j]<<' '<<chromosome[i].y[j]<<' ';
    }
    cout<<endl;
}
cout<<endl;
*/
     int    ticks = 1;
     double value = des[0].coverage;   //ticks为保留直接遗传的染色体
     for(int i=1;i<chromosome_num / 2;i++){
            if(des[i].coverage != value) break;
            ticks++;
     }
     double coverage_t[chromosome_num];   //保存以便还原
     for(int i=0;i<chromosome_num;i++) coverage_t[i] = des[i].coverage;
     for(int i=0;i<chromosome_num;i++) des[i].coverage -= (des[chromosome_num - 1].coverage - 0.01);
     //最后加上0.01是防止des[i].coverage全变为0，sum为0导致除法出错
     double sum = 0, select_p, accu[chromosome_num + 1];   //accu为累积选中概率

     for(int i=ticks;i<chromosome_num;i++) {
            sum += des[i].coverage;
            accu[i] = 0;
     }

     accu[chromosome_num] = 0;
     Chromosomes temp[chromosome_num];
     int count_num = 0;
     for(;count_num < ticks;count_num++)
        for(int c=0;c<chromosome[count_num].x.size();c++){
            temp[count_num].x.push_back(chromosome[count_num].x[c]);
            temp[count_num].y.push_back(chromosome[count_num].y[c]);
        }

     for(int i=ticks + 1;i<=chromosome_num;i++){
        select_p = des[i-1].coverage / sum;
        accu[i]  = accu[i-1] + select_p;
//cout<<des[i-1].coverage<<' '<<select_p<<' '<<accu[i]<<endl;
     }

//for(int i=0;i<=chromosome_num;i++) cout<<accu[i]<<' ';cout<<endl;
     for(int i=ticks;i<chromosome_num;i++){     //选chromosome_num - ticks次
        double p = (double)(rand() * 1.0 /(RAND_MAX+1.0));
        for(int j=ticks;j<chromosome_num;j++){
            if(p >= accu[j] && p <= accu[j+1]){
//cout<<j<<' ';
//cout<<"1111111111111111111111"<<endl;
                for(int c=0;c<chromosome[j].x.size();c++){
                     temp[count_num].x.push_back(chromosome[j].x[c]);
                     temp[count_num].y.push_back(chromosome[j].y[c]);
                }
//cout<<chromosome[j].x.size()<<' '<<chromosome[j].y.size()<<' '<<temp[count_num].x.size()<<' '<<temp[count_num].y.size()<<endl;
                count_num++;
                break;
            }
        }
//cout<<des[i-1].coverage<<' '<<select_p<<' '<<accu[i]<<endl;
     }
//cout<<endl;
//cout<<count_num<<endl;
        //选择概率，等于适应值大小除以适应值总和
     for(int i=0;i<chromosome_num;i++){
        for(int j=0;j<chromosome[i].x.size();j++){
            chromosome[i].x[j] = temp[i].x[j];
            chromosome[i].y[j] = temp[i].y[j];
//cout<<chromosome[i].x[j]<<' '<<chromosome[i].y[j]<<' ';
        }
//cout<<endl;
     }
//cout<<ticks<<' ';
//cout<<endl;cout<<endl;
     //sort(des, des + chromosome_num, cmp2);
     for(int i=0;i<chromosome_num;i++) des[i].coverage = coverage_t[i];
}

//评估染色体的适应度
void evaluate()
{
     Node node_t;
     for(int i=0;i<chromosome_num;i++){
        for(int i=0;i<grid_num;i++) cover_flag[i] = static_cover[i];
//for(int i=0;i<grid_num;i++) cout<<cover_flag[i]<<' ';cout<<endl;
        //每次评估前均作复原
         for(int j=0;j<chromosome[i].x.size();j++){
              node_t.x = chromosome[i].x[j];
              node_t.y = chromosome[i].y[j];
//cout<<node_t.x<<' '<<node_t.y<<' ';
              cover(node_t);
         }
//cout<<endl;
//for(int i=0;i<grid_num;i++) cout<<cover_flag[i]<<' ';cout<<endl;
         des[i].coverage = compute_coverage();
         if(max_coverage < des[i].coverage){
                max_coverage = des[i].coverage;
                move_d       = des[i].move_d;
         }
         des[i].init_c = des[i].coverage;

         double move_dis = 0.0;
         int mobile_num = chromosome[0].x.size();
         for(int k=0;k<mobile_num;k++){
              int dx = node[mobile_pos[k]].x - chromosome[i].x[k];
              int dy = node[mobile_pos[k]].y - chromosome[i].y[k];
              move_dis += sqrt(dx * dx + dy * dy);
         }
         move_dis /= mobile_num;
         des[i].move_d = move_dis;
/*
//cout<<des[i].coverage<<' ';
         des[i].coverage = des[i].coverage * 1000 / (move_dis + 100);
         //覆盖率综合移动距离作为后代选择标准，慎选
//cout<<des[i].coverage<<endl;
*/

//cout<<move_dis<<endl;
         des[i].label = i;
     }
/*
    //下面做sigma变换（标准化），目的是为了方便计算选择概率
    double sum = 0.0, average, deviation;   //deviation为样本方差
    for(int i=0;i<chromosome_num;i++)  sum += des[i].coverage;
    average = sum / chromosome_num;     //适应值的平均数
    sum = 0.0;
    for(int i=0;i<chromosome_num;i++)  sum += pow(des[i].coverage - average, 2.0);
    deviation = sum / (chromosome_num - 1) * 1.0;
    for(int i=0;i<chromosome_num;i++) des[i].coverage = 1.0 + (des[i].coverage - average) / 2.0 * deviation;
*/
    sort(des, des + chromosome_num, cmp1);
    //对des以覆盖率排序
    for(int i=0;i<chromosome_num;i++){
        int t = des[i].label;
        chromosome[t].label = i;
    }
}//evaluate

double evaluate_single(Chromosomes chro)
{
     //extern Chromosomes chromosome[DEF_CHROMOSOME_NUM];        //染色体
     Node node_t;
     for(int i=0;i<grid_num;i++) cover_flag[i] = static_cover[i];
     for(int j=0;j<chro.x.size();j++){
            node_t.x = chro.x[j];
            node_t.y = chro.y[j];
            cover(node_t);
     }
     return compute_coverage();
}

//随机生成染色体
void produce_chrom()
{
//std::cout<<field[0]<<' ' <<field[2]<<std::endl;
     srand((unsigned)time(NULL));
     for(int i=0;i<chromosome_num;i++){
         for(int j=0;j<mobile_num;j++){
            chromosome[i].x.push_back(rand()%(field[1] - field[0] + 1) + field[0]);
            chromosome[i].y.push_back(rand()%(field[3] - field[2] + 1) + field[2]);
         }
     }
}//produce_chrom

//计算覆盖率
double compute_coverage()
{
    int count_nums = 0;
    for(int i=0;i<grid_num;i++)
        if(cover_flag[i])  count_nums++;
    return count_nums * 1.0 / grid_num;
}

//随机生成移动节点并计算静态节点的覆盖率
void static_covere()
{
    for(int i=0;i<grid_num;i++)  cover_flag[i] = false;
    for(int i=0;i<sensor_num;i++)  is_mobile[i] = false;
    for(int i=0;i<mobile_num;i++){         //生成移动节点
        int a = (i-1) * sensor_num / mobile_num + 1;
        int b = i * sensor_num / mobile_num + 1;
        int k;
        do{
           k = rand()%(b - a + 1) + a;
        }while(is_mobile[k]);
        is_mobile[k] = true;            //结点k为移动节点
        mobile_pos[i] = k;
    }
    for(int i=0;i<sensor_num;i++)
        if(!is_mobile[i]) cover(node[i]);
    printf("static sensor coverage: %.2f%%\n",100.0 * compute_coverage());

}

//对传感器节点node进行范围覆盖
void cover(Node node2)
{
    double minx = node2.x - sense_range;    //感知半径范围x
    minx = minx < field[0] ? field[0] : minx;
    double maxx = node2.x + sense_range;
    maxx = maxx > field[1] ? field[1] : maxx;
    double miny = node2.y - sense_range;    //感知半径范围y
    miny = miny < field[2] ? field[2] : miny;
    double maxy = node2.y + sense_range;
    maxy = maxy > field[3] ? field[3] : maxy;

    double d = (minx - field[0]) * 1.0 / grid;
    int x1 = d - (int)d > 0.5 ? (int)d + 1 : (int)d;    //左端离左边界的格数
    d = (maxx - field[0]) * 1.0 / grid;
    int x2 = d - (int)d > 0.5 ? (int)d + 1 : (int)d;    //右端离左边界的格数
    d = (miny - field[2]) * 1.0 / grid;
    int y1 = d - (int)d > 0.5 ? (int)d + 1 : (int)d;    //上端离上边界的格数
    d = (maxy - field[2]) * 1.0 / grid;
    int y2 = d - (int)d > 0.5 ? (int)d + 1 : (int)d;    //下端离上边界的格数

    double x, y, dis, des = (double)pow(sense_range,2.0);
    for(int i=x1;i<x2;i++){
        for(int j=y1;j<y2;j++){
            x   = field[0] + (i + 0.5) * grid * 1.0;
            y   = field[2] + (j + 0.5) * grid * 1.0;
            dis = (double)pow(abs(x - node2.x),2.0) + (double)pow(abs(y - node2.y),2.0);
            //printf("%f %f\n",dis,des);
            if(dis < des)
                cover_flag[col_grid_num * j + i] = true;
        }
    }
}

void cover_first()
{
//std::cout<<sensor_num<<std::endl;
    for(int i=0;i<sensor_num;i++)  cover(node[i]);
}

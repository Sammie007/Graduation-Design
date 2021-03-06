#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QDebug>
#include <QTime>
#include <QDir>
#include <QFile>
#include <iostream>

//#pragma execution_character_set("utf-8")

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{

}

void MainWindow::on_radioButton_clicked()   //随机产生
{

}

void MainWindow::on_radioButton_2_clicked()     //文件导入
{
    extern int read_file(const char * const filename);
    extern int sensor_num;
    QString filename = QFileDialog::getOpenFileName(this, tr("选择文件"), "./", tr("Text files(*.txt)"));
    //std::cout<<filename.toStdString()<<std::endl;
    char *file_name = (char*)filename.toStdString().c_str();
    sensor_num = read_file(file_name);
//std::cout<<"finish"<<std::endl;
}

void MainWindow::on_pushButton_2_clicked()
{
    dialog = new Dialog(this);
    dialog->setModal(false);
    dialog->show();

    clock_t start,finish;      //程序计时
    double totaltime = 0.0;
    start = clock();

    extern int  field[4];
    field[0] = 190;  field[1] = field[0] + 280;
    field[2] = 84;   field[3] = field[2] + 280;

    extern int     grid;
    extern int     row_grid_num; row_grid_num = (field[3]-field[2]) / grid;
    extern int     col_grid_num; col_grid_num = (field[1]-field[0]) / grid;
    extern int     grid_num;     grid_num = row_grid_num * col_grid_num;     //平面上单元格的数目
//std::cout<<grid_num<<std::endl;
    extern int     placement_unit; placement_unit = (field[1] - field[0]) / 10;
    extern int     max_move_dist; max_move_dist = 5 * placement_unit;
    extern double  borderlength;  borderlength = 0.8 * (field[1]-field[0]);
    extern double  mobile_percent; extern int  sensor_num;
    extern int     mobile_num;   mobile_num= mobile_percent * sensor_num;
    extern double  transmission_range; transmission_range = 2.0 * placement_unit * sqrt(60.0 / sensor_num);
    extern double  sense_range;   sense_range = transmission_range / 2.0;
//std::cout<<sensor_num<<std::endl;
    extern void produce_chrom();
    produce_chrom();

    extern bool  cover_flag[];       //两个状态，1-覆盖，0-未覆盖
    for(int i=0;i<grid_num;i++)  cover_flag[i] = false;  //初始化
    //extern void cover();
    extern void cover_first();
    cover_first();
    //对所有传感器节点进行初始覆盖

    extern double compute_coverage();
    double max_goal = compute_coverage();    //计算初始覆盖率
    printf("initial coverage: %.2f%%\n",100.0 * max_goal);      //并输出

    extern void static_covere(); static_covere();
    //随机生成移动节点并计算静态节点的覆盖率
    extern bool static_cover[];
    for(int i=0;i<grid_num;i++) static_cover[i] = cover_flag[i];    //保存静态节点的覆盖
    //extern chromo_attribute des[];
    extern double ideal_length; ideal_length = sqrt((field[1]-field[0])*(field[3]-field[2]) / sensor_num);
    //理想边长
    extern double max_coverage; max_coverage = 0.0;
    extern double move_d; move_d = 0.0;    //记录最大覆盖率和移动距离

    extern int max_generations;
    extern void evaluate();
    extern void select2();
    extern void single_point_crossover();
    extern void variation(int generation);
    extern int  chromosome_num;
    extern double probability_crossover;
    extern double probability_mutate;
    //extern chromo_attribute des[];
    for(int generation = 0; generation < max_generations; generation++){       //遗传算法主体
         evaluate();
         if(max_coverage == 1.0) break;
//for(int i=0;i<chromosome_num;i++) printf("%.2f ",des[i].init_c); printf("\n");
         select2();      //选择机制
         //选择决定了个体的优胜劣汰
         single_point_crossover();
         //交叉保证了优质个体的基因重组
         variation(generation);
         //变异是保证不陷入局部最优的重要保证
    }

    //printf("选择策略： 轮盘赌选法\n");
    //printf("杂交策略： 算术单点交叉\t");
    printf("probability of crossover: %d%%\n",(int)(probability_crossover * 100));
    //printf("变异策略： 单点邻域变异\t");
    printf("probability of mutate: %d%%\n",(int)(probability_mutate * 100));
    printf("\nCount:\nconvergence generations: %d\n",max_generations);
    printf("\nchromosome nums: %d\n",chromosome_num);
    printf("\nmax_coverage: %f%%\n",max_coverage * 100);
    printf("\nmove_distance: %.2f\n",move_d);
    //for(int i=0;i<chromosome_num;i++) printf("%d\n",des[i].is_select);
    finish = clock();
    totaltime = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("\ncost times: %f Seconds\n",totaltime);
}

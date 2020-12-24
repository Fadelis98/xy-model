#include "MyLattice.h"
#include "anyMc.h"
#include <fstream>
#include <iostream>

using namespace lat;
using namespace std;
namespace lat
{
    static std::random_device __xymodel_rd;
    static std::default_random_engine __xychoose_e(__xymodel_rd());
    template <>
    void chooseRvec(double &numOfpi)
    {
        static uniform_real_distribution<double> __xychoose_ur(0, 2);
        numOfpi = __xychoose_ur(__xychoose_e);
    }

    template <>
    double reflect(const double &orign, const double &rvec)
    {
        double val = rvec - (orign - rvec);
        if (val >= 2)
        {
            return val - 2;
        }
        else if (val < 0)
        {
            return val + 2;
        }
        return val;
    }

    template <>
    double caculatePairEnergy(const double &sigma1, const double &sigma2)
    {
        return -cos((sigma1 - sigma2)*M_PI);
    }

    template <>
    double caculateMagnitude(const Lattice2D<double> &model)
    {
        int scale = model.Scale();
        double MagX, MagY = 0;
        for (size_t x = 0; x < scale; x++)
        {
            for (size_t y = 0; y < scale; y++)
            {
                double deg = model.Get(x, y) * M_PI;
                MagX += cos(deg);
                MagY += sin(deg);
            }
        }
        return sqrt(MagX * MagX + MagY * MagY);
    }
} // namespace lat
void outputLattice(const Lattice2D<double> &model,string filename);
void outputVortices(const vector<vector<position2d>> &vorts,string filename);
double relevantfunc(const Lattice2D<double>& model,int distance,double beta);
vector<vector<position2d>> searchVortices(const Lattice2D<double>& model,int r);
double existVortices(const Lattice2D<double> &model,const position2d& p,int r);
vector<position2d> loopGene(const position2d& center,int r);


inline void showVortsOnLat(XY2D model,double temperature,int scale,vector<int> VscaleList)
{
    outputLattice(model,"LatOf"+to_string(temperature)+"s"+to_string(scale)+".csv");
    for (auto v:VscaleList)
    {
        auto vorts =  searchVortices(model,v);
        for(auto vort:vorts)
        {
            if (vort.size()!=0)
            {
                outputVortices(vorts,"VortOf"+to_string(temperature)+"s"+to_string(scale)+ "v"+to_string(v) + ".csv");
                break;
            }
        }         
    }
}



vector<double> TList = {0.89,0.9};
vector<int> ScaleList = {64};
vector<int> VortScaleList = {1,2,3,5,8,16,24 };

int main()
{
    int numOfrel = 0;
    int numOfdataPoints = 0;
    int preheat_iter = 3000;
    int intervalWolff = 100;
    int intervalMP = 1000;
    for (auto Latscale : ScaleList)
    {
        for (auto temperature : TList)
        {
            double beta = 1 / temperature;
            XY2D wolffmodel(Latscale);
            for (size_t i = 0; i < preheat_iter; i++) //预热模型
            {
                WolffProcess(wolffmodel, beta);
            }
            for (size_t i = 0; i < preheat_iter; i++) //预热模型
            {
                for (size_t j = 0; j < 100; j++)
                {
                    MetroPolisProcess(wolffmodel, beta);
                } 
            }
            showVortsOnLat(wolffmodel,temperature,Latscale,VortScaleList);

        }

    }
}


/*Mag
vector<double> TList = linspace(0,3,31);
vector<int> ScaleList = {20};
int main()
{
    int numOfdataPoints = 2000;
    int preheat_iter = 0;
    int intervalWolff = 10;
    int intervalMP = 100;
    for (auto scale : ScaleList)
    {
        for (auto temperature : TList)
        {
            double beta = 1 / temperature;

            XY2D MPmodel(scale);
            ofstream dataofMP("MPmag" + to_string(temperature) + ".csv");
            for (size_t i = 0; i < preheat_iter; i++) //预热模型
            {
                MetroPolisProcess(MPmodel, beta);
            }
            for (size_t i = 0; i < numOfdataPoints; i++) //采样
            {
                for (size_t j = 0; j < intervalMP; j++) //采样间隔
                {
                    MetroPolisProcess(MPmodel, beta);
                }
                double mag = caculateMagnitude(MPmodel);
                dataofMP << mag << endl;
            }
            dataofMP.close();

            XY2D Wolffmodel(scale);
            ofstream dataofWolff("Wolffmag" + to_string(temperature) + ".csv");
            for (size_t i = 0; i < preheat_iter; i++) //预热模型
            {
                MetroPolisProcess(Wolffmodel, beta);
            }
            for (size_t i = 0; i < numOfdataPoints; i++) //采样
            {
                for (size_t j = 0; j < intervalWolff; j++) //采样间隔
                {
                    MetroPolisProcess(Wolffmodel, beta);
                }
                double mag = caculateMagnitude(Wolffmodel);
                dataofWolff << mag << endl;
            }
            dataofMP.close();

        }
    }
}
*/

/*
vector<double> TList = linspace(0, 6, 31);
vector<int> ScaleList = {64};
int main()
{
    int numOfrel = 10;
    int numOfdataPoints = 1000;
    int preheat_iter = 3000;
    int intervalWolff = 200;
    int intervalMP = 1000;
    for (auto scale : ScaleList)
    {

        ofstream dataOfRel("Rel" + to_string(scale) + ".csv");
        dataOfRel << "temp";
        for (size_t i = 1; i <= numOfrel; i++)
        {
            dataOfRel << ",d=" << i;
        }
        dataOfRel << endl;
        for (auto temperature : TList)
        {
            double beta = 1 / temperature;
            #pragma region 定义要计算的物理量相关变量
            ofstream dataofWolff("Wolffmag" + to_string(temperature) + ".csv");
            vector<vector<double>> Rel(numOfrel,vector<double>(numOfdataPoints));
            #pragma endregion

            #pragma region Wolff算法相关段落
            XY2D wolffmodel(scale);
            for (size_t i = 0; i < preheat_iter; i++) //预热模型
            {
                WolffProcess(wolffmodel, beta);
            }
            for (size_t i = 0; i < numOfdataPoints; i++) //采样
            {
                for (size_t j = 0; j < intervalWolff; j++) //采样间隔
                {
                    WolffProcess(wolffmodel, beta);
                }
                #pragma region 计算模型相关物理量
                
                for (size_t j = 1; j <= numOfrel; j++)//计算不同距离的关联函数
                {
                    Rel[j-1][i] = relevantfunc(wolffmodel,j,beta);
                }
                #pragma endregion 计算模型相关物理量
            }
            #pragma endregion wolff算法
            dataOfRel << temperature;
            for (auto reld:Rel)
            {
                dataOfRel  << "," << mean(reld);
            }
            dataOfRel << endl;
        }
        dataOfRel.close();
    }
}
*/
#pragma region 计算关联函数
double relevantfunc(const Lattice2D<double>& model,int distance,double beta)
{
    int scale = model.Scale();
    vector<double> rel;
    rel.reserve(scale*scale*2);
    for (size_t i = 0; i < scale; i++)
    {
        for (size_t j = 0; j < scale; j++)
        {
            rel.push_back(lat::caculatePairEnergy(model.Get(i,j),model.Get(i+distance,j)));
            rel.push_back(lat::caculatePairEnergy(model.Get(i,j),model.Get(i,j+distance)));
        }
    }
    return -lat::mean(rel);
}
#pragma endregion

#pragma region 寻找拓扑激发
vector<vector<position2d>> searchVortices(const Lattice2D<double>& model,int r)
{
    int scale = model.Scale();
    vector<position2d> vortices;
    vector<position2d> antiVortices;
    vector<position2d> doubleVortices;
    vector<position2d> doubleAntivortices;
    int interval = 1;
    for (int i = 0; i < scale; i+=interval)
    {
        for (int j = 0; j < scale; j += interval)
        {
            double typeOfP = existVortices(model,position2d{i,j},r);
            if(typeOfP >= 0.999)
            {
                vortices.push_back(position2d{i,j});
            }
            else if (typeOfP >= 1.999)
            {
                doubleVortices.push_back(position2d{i,j});
            }
            else if (typeOfP <= -0.999)
            {
                antiVortices.push_back(position2d{i,j});
            }
            else if (typeOfP <= -1.999)
            {
                doubleAntivortices.push_back(position2d{i,j});
            }
          
        }
    }
    return vector<vector<position2d>>{vortices,antiVortices,doubleVortices,doubleAntivortices};

}

double existVortices(const Lattice2D<double> &model,const position2d& p,int r)
{
    auto loop =  loopGene(p,r);
    position2d lastNode;
    position2d newNode= loop[0];
    double delta=0;
    for(auto loopNode:loop)
    {
        lastNode = newNode;
        newNode = loopNode;
        double d = model.Get(newNode.x,newNode.y)-model.Get(lastNode.x,lastNode.y);
        if (d>1)
        {
            d -= 2;

        }
        else if (d<-1)
        {
            d += 2;
        }
        if(abs(d) >0.45)
        {
            return 0;
        }
        delta += d;
    }
    return delta/2;
}

vector<position2d> loopGene(const position2d& center,int r)
{
    vector<position2d> loop;
    loop.reserve(8*r);
    for (int i = -r+1; i <= r; i++)
    {
        loop.push_back(position2d{center.x+i,center.y-r});
    }
    for (int i = -r+1; i <= r; i++)
    {
        loop.push_back(position2d{center.x+r,center.y+i});
    }
        for (int i = -r+1; i <= r; i++)
    {
        loop.push_back(position2d{center.x-i,center.y+r});
    }
        for (int i = -r+1; i <= r; i++)
    {
        loop.push_back(position2d{center.x-r,center.y-i});
    }
    loop.push_back(position2d{center.x+1-r,center.y-r});
    return loop;
}


#pragma endregion

#pragma region 输出数据
void outputLattice(const Lattice2D<double> &model,string filename)
{
    std::ofstream dataOfLat(filename);
    for(auto& line:model.Data())
    {
        for(auto &p:line)
        {
            dataOfLat << p << "," ;
        }
        dataOfLat << endl;
    }
    dataOfLat.close();
}

void outputVortices(const vector<vector<position2d>> &vorts,string filename)
{
    std::ofstream dataOfVort(filename);
    for(auto& line:vorts)
    {
        for(auto &p:line)
        {
            dataOfVort  <<p.x << "," <<p.y  << "," ;
        }
        dataOfVort << endl;
    }
    dataOfVort.close();
}

#pragma endregion
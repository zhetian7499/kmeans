#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <omp.h>
using namespace std;
#define DIMENSION 2
#define KCLUSTER 3
#define N 20
#define ITERATOR 50

class Point{
private:
    int numCluster;
    vector<double> values;

public:
    Point(vector<double> value, int init = 0):values(value), numCluster(init){}

    int get_numCluster() { return numCluster; }
    //something needed improve
    bool equal(Point y) { return values == y.values; }

    void set_numCluster(int val) { numCluster = val; }

    vector<double> getVec() { return values; }
    
    double calcDistance(vector<double> other){
        double dist;
        #pragma omp parallel for reduction(+: dist) num_threads(16)
        for (int i = 0; i < other.size(); i++)
		    dist += pow(values[i] - other[i], 2);
	    return sqrt(dist);
    }
};

class Cluster{
private:
    vector<double> centre;
    vector<Point> points;

public:
    Cluster(int i, Point centre){
        this->centre = centre.getVec();
        this->addPoint(centre);
        centre.set_numCluster(i);
    }

    Point getPoint(int pos) { return points[pos]; }

    int getSize() { return points.size(); }

    void addPoint(Point p){
        points.push_back(p);
    }

    void removePoint(Point point){
        int size = points.size();
        #pragma omp parallel for num_threads(16)
        for (int i = 0; i < size; i++){
            if (points[i].equal(point))
                points.erase(points.begin() + i);
        }
    }
    
    void removeAllpoints() { points.clear(); }

    vector<double> getCentre() { return centre; }

    void setCentre(int pos, double val) { this->centre[pos] = val; }
};

class KMeans{
private:
    int k, iters;
    vector<Cluster> clusters;

    int getNearestCluster(Point point, vector<Cluster> clusters){
        double min_dist;
        int NearestCluster;

        min_dist = point.calcDistance(clusters[0].getCentre());
        NearestCluster = 1;

        #pragma omp parallel for num_threads(16)
        for (int i = 1; i < k; i++){
            double dist;
            dist = point.calcDistance(clusters[i].getCentre());
            if (dist < min_dist){
                min_dist = dist;
                NearestCluster = i + 1;
                
            }
            // cout << i << " ? " << min_dist << "  " << dist << endl;
            // cout << "pass" << endl;
        }
        // cout << "this ended" << endl;
        return NearestCluster;
    }

public:
    KMeans(int k_number, vector<Point> &all_points, int iter = ITERATOR): k(k_number), iters(iter){
        check(all_points);
        run(all_points); 
        output_Points(all_points);
        output_Clusters();
    }

    //ouput warnings
    void check(vector<Point> &all_points){
        //if number of clusters > number of points
        if(all_points.size() < KCLUSTER){
            cout << "Number of clusters is smalled than number of points." << endl;
        }
        //if dimension < 2
        if(DIMENSION < 2){
            cout << "Number of DIMENSION is below 2." << endl; 
        }
        //if number of cluster < 2
        if(DIMENSION < 2){
            cout << "Number of cluster is below 2." << endl; 
        }
    }
    
    void initCluster(vector<Point> &all_points){
        vector<int> used_pointIds;
        for (int i = 1; i <= k; i++){
            while (1){
                int index = rand() % N;
                if (find(used_pointIds.begin(), used_pointIds.end(), index) == used_pointIds.end()){
                    used_pointIds.push_back(index);
                    all_points[index].set_numCluster(i);
                    Cluster cluster(i, all_points[index]);
                    clusters.push_back(cluster);
                    break;
                }
            }
        }
        //cout << "init successfully" << endl
    }

    //used to relocate centre
    void recalCentre(){
        for (int i = 0; i < k; i++){
            int ClusterSize = clusters[i].getSize();

            for (int j = 0; j < DIMENSION; j++){
                double sum = 0.0;
                if (ClusterSize > 0){
                    #pragma omp parallel for reduction(+: sum) num_threads(16)
                    for (int p = 0; p < ClusterSize; p++){
                        sum += clusters[i].getPoint(p).getVec()[j];
                    }
                    clusters[i].setCentre(j, sum / ClusterSize);
                }
            }
        }
    }

    void output_Points(vector<Point> &all_points){
        ofstream pointsFile;
        pointsFile.open("points.csv", ios::out);

        for (int i = 0; i < N; i++){
            pointsFile << all_points[i].getVec()[0];
            for (int j = 1; j < DIMENSION; j++){
                pointsFile << ',' << all_points[i].getVec()[j];
            }
            pointsFile << ',' <<all_points[i].get_numCluster() << endl;
        }
        pointsFile.close();
    }

    void output_Clusters(){
        ofstream clusterFile;
        clusterFile.open("clusters.csv");

        for (int i = 0; i < k; i++){
            //cout << "Cluster " << i + 1 << " centre : ";
            //cout << clusters[i].getCentre()[0] << " ";
            clusterFile << clusters[i].getCentre()[0];
            for (int j = 1; j < DIMENSION; j++){
                //cout << clusters[i].getCentre()[j] << " ";    
                clusterFile << ',' << clusters[i].getCentre()[j];
            }
            //cout << endl;
            clusterFile << endl;
        }
        clusterFile.close();
    }
    
    void run(vector<Point> &all_points){
        initCluster(all_points);
        
        int iter = 1;
        while (true){
            //cout << iter << " times" << endl;
            bool done = true;

            // resign the rest points
            #pragma omp parallel for num_threads(16)
            for (int i = 0; i < N; i++){
                int currentCluster = all_points[i].get_numCluster();
                int nearestCluster = getNearestCluster(all_points[i], clusters);
                if (currentCluster != nearestCluster){
                    all_points[i].set_numCluster(nearestCluster);           
                    done = false;
                }
            }
           
            //all clusters are ruined
            #pragma omp parallel for num_threads(KCLUSTER)
            for (int i = 0; i < k; i++)
                clusters[i].removeAllpoints();

            // reassign all_points
            for (int i = 0; i < N; i++){
                // cluster index is ID-1
                clusters[all_points[i].get_numCluster() - 1].addPoint(all_points[i]);
            }
            
            recalCentre();            

            if (done || iter >= iters){
                //cout << "completed by " << iter << " times" << endl;
                cout << "Please check data in .csv in this folder." << endl;
                break;
            }
            iter++;
        }
    }
};
    
class Data{
public:
    Data(int n){
        points = initPoints();
    }
    vector<Point> getPoints(){ return points; }

private:
    vector<Point> points;  
    vector<Point> initPoints(){
        vector<Point> all_points;

        random_device rd;  //obtain a seed for the random number engine
	    mt19937 gen(rd()); 
	    std::uniform_real_distribution<float> distribution(0.0, 30.0);
    
        for (int i = 0; i < N; i++){
            vector<double> value;
            for (int j = 0 ; j < DIMENSION; j++){
                value.push_back(distribution(gen));
            }
            Point point(value);
            all_points.push_back(point);
        }
        return all_points;
    } 
};

int main(){
    //prepare random data
    vector<Point> all_points = Data(N).getPoints();
    if (all_points.size() == N)
        cout << "Data are well prepared."<<endl;
    
    //run K-means
    KMeans kmeans(KCLUSTER, all_points);
    return 0;
}

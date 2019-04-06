#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <omp.h>

#define CHUNK 1
#define THREAD_NUM 64

class Coordinate
{
    public:
	float x;
	float y;
	float z;

	Coordinate operator+(const Coordinate& coor){
		Coordinate newcoor;
		newcoor.x = this->x + coor.x;
		newcoor.y = this->y + coor.y;
		newcoor.z = this->z + coor.z;
		return newcoor;
	}
	Coordinate operator-(const Coordinate& coor){
		Coordinate newcoor;
		newcoor.x = this->x - coor.x;
		newcoor.y = this->y - coor.y;
		newcoor.z = this->z - coor.z;
		return newcoor;
	}
	double operator*(const Coordinate& coor){
		return (this->x * coor.x + this->y * coor.y + this->z * coor.z);
	}
};

Coordinate Cross (const Coordinate& coor1, const Coordinate& coor2);
double Size(const Coordinate& coor1);

void PrintUsage();

int main (int argc, char *argv[])
{
	char *inputFile = argv[1];
	char *outputFile = argv[2];

	//int interval = atoi(argv[3]);
	//char *para_1;
        //int  interval = strtol(argv[3], &para_1, 10);
        int interval = 1;

        if(argc != 3){
                PrintUsage();
                return EXIT_FAILURE;
        }

	std::ifstream Fin(inputFile);
	if(!Fin) std::cerr<<"Cannot open input file " << inputFile << std::endl;

	std::vector<Coordinate> atoms;
	std::string buf;
	Coordinate xyz;
	
	int number = 1;
	while (getline(Fin, buf)){
		if(number % interval == 0){
			xyz.x = atof(buf.substr(31, 11).c_str());
			xyz.y = atof(buf.substr(43, 11).c_str());
			xyz.z = atof(buf.substr(55, 11).c_str());
			atoms.push_back(xyz);
		}
		number++;
	}
	atoms.push_back(atoms[0]);

	Fin.close();
	std::cout << "Total atom number: " << atoms.size()-1 << std::endl;

	double W1 = 0;
	double size = 0;
	double W2 = 0;
	double total = 0;
	double ptotal = 0;

    	int i, j;
    	Coordinate s1, s2, V, C;

    	#pragma omp parallel private(j, s1, s2, V, C, W1, size, W2, ptotal) reduction(+:total) 
    	{
		#pragma omp for schedule(dynamic, CHUNK)
        	for(i=0; i<atoms.size()-1; i++) {
			for(j=i+1; j<atoms.size()-1; j++) {

				s1 = atoms[i+1] - atoms[i];		
				s2 = atoms[j+1] - atoms[j];		

				V = atoms[i] -  atoms[j];
				C = Cross(s1, s2);

				W1 = C * V;
				size = Size(V);
				W2 = W1 / (size * size * size);
				ptotal += W2;
			}
			if(i % 10000 == 0) {
				std::cout << inputFile << " finish:" << i << std::endl;
			}
        	}
		total += ptotal;
    	}

	total = total / (2 * M_PI); 

	std::ofstream Fout(outputFile);
	if(!Fout)std::cerr<<"Cannot open output file " << outputFile << std::endl;

	
        std::cout << inputFile << " Calculated writhe number: " 
	<< std::setprecision(3) << std::fixed
	<< total <<  "\n";

        Fout << inputFile << " Calculated writhe number: " 
	<< std::setprecision(3) << std::fixed
	<< total <<  "\n";
        Fout.close();
        return 0;
}

Coordinate Cross (const Coordinate& coor1, const Coordinate& coor2)
{
	Coordinate newcoor;
	newcoor.x = coor1.y * coor2.z - coor1.z * coor2.y;
	newcoor.y = coor1.z * coor2.x - coor1.x * coor2.z;
	newcoor.z = coor1.x * coor2.y - coor1.y * coor2.x;
	return newcoor;
}

double Size(const Coordinate& coor1)
{
	return sqrt(coor1.x*coor1.x + coor1.y*coor1.y + coor1.z*coor1.z);
}

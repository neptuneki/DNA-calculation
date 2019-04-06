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

#define ZERO 1.0e-5

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

	Coordinate operator/(const double &a){
		Coordinate newcoor;
		newcoor.x = this->x / a;
		newcoor.y = this->y / a;
		newcoor.z = this->z / a;
		return newcoor;
	}

	Coordinate operator*(const double &a){
		Coordinate newcoor;
		newcoor.x = this->x * a;
		newcoor.y = this->y * a;
		newcoor.z = this->z * a;
		return newcoor;
	}

	double operator*(const Coordinate& coor){
		return (this->x * coor.x + this->y * coor.y + this->z * coor.z);
	}

};

const Coordinate Cross (const Coordinate& coor1, const Coordinate& coor2);

const double Link (Coordinate& t1, Coordinate& t2, Coordinate& rjoin);

const double Size(const Coordinate& coor1);

void PrintUsage();

int main (int argc, char *argv[])
{
	char *inputFile = argv[1];
	char *outputFile = argv[2];

        if(argc != 3){
                PrintUsage();
                return EXIT_FAILURE;
        }

	std::ifstream Fin(inputFile);
	if(!Fin) std::cerr<<"Cannot open input file " << inputFile << std::endl;

	std::vector<Coordinate> PDB1;
	std::vector<Coordinate> PDB2;
	Coordinate xyz;
	std::string buf;
	
	int number = 1;

	while (getline(Fin, buf)){
	    xyz.x = atof(buf.substr(31, 11).c_str());
	    xyz.y = atof(buf.substr(43, 11).c_str());
	    xyz.z = atof(buf.substr(55, 11).c_str());

	    if(number%2 == 1)  {
		PDB2.push_back(xyz);
	    }
	    else  {
		xyz.x = 0.5 * (xyz.x + PDB2.back().x);
		xyz.y = 0.5 * (xyz.y + PDB2.back().y);
		xyz.z = 0.5 * (xyz.z + PDB2.back().z);
		PDB1.push_back(xyz);
	    }
	    number++;
	}

	PDB1.push_back(PDB1[0]);
	PDB2.push_back(PDB2[0]);

	Fin.close();
	std::cout << "Total bps: " << PDB1.size()-1 << " should be equal to " << PDB2.size()-1 << std::endl;

	double total = 0;
	double ptotal = 0;
    	int i, j;
	double S1, S2;
	Coordinate segment1, segment2, r, rtemp, k1, k2;

    	#pragma omp parallel private(j, ptotal, S1, S2, segment1, segment2, r, rtemp, k1, k2) reduction(+:total) 
    	{
		#pragma omp for schedule(dynamic, CHUNK)
        	for(i=0; i<PDB1.size()-1; i++) {
		    for(j=0; j<PDB2.size()-1; j++) {
				
        		segment1 = PDB1[i+1] - PDB1[i];
        		segment2 = PDB2[j+1] - PDB2[j];

			S1 = Size(segment1);
			S2 = Size(segment2);

        		k1 = segment1 / S1;
        		k2 = segment2 / S2;

			rtemp = PDB1[i] - PDB2[j];


			r = rtemp + segment1;
			ptotal += Link(k1, k2, r);

			r = rtemp - segment2;
			ptotal += Link(k1, k2, r);

			r = rtemp + segment1 - segment2;
			ptotal -= Link(k1, k2, r);

			r = rtemp;
			ptotal -= Link(k1, k2, r);

		    }
		    if(i % 10000 == 0) {
			std::cout << inputFile << " finish:" << i << std::endl;
		    }
        	}
		total += ptotal;
    	}

	total = total / (4 * M_PI); 

	std::ofstream Fout(outputFile);
	if(!Fout)std::cerr<<"Cannot open output file " << outputFile << std::endl;
	
        std::cout << std::setprecision(3) << std::fixed << inputFile << " Calculated linking number: " << total <<  "\n";
        Fout << std::setprecision(3) << std::fixed << inputFile << " Calculated linking number: " << total <<  "\n";
        Fout.close();
        return 0;
}

void PrintUsage()
{
        std::cout << "\nUsage: XXX.exe <input_file> <output_file>\n\n";

}

const Coordinate Cross (const Coordinate& coor1, const Coordinate& coor2)
{
	Coordinate newcoor;
	newcoor.x = coor1.y * coor2.z - coor1.z * coor2.y;
	newcoor.y = coor1.z * coor2.x - coor1.x * coor2.z;
	newcoor.z = coor1.x * coor2.y - coor1.y * coor2.x;
	return newcoor;
}

const double Size(const Coordinate& coor1)
{
	return sqrt(coor1.x*coor1.x + coor1.y*coor1.y + coor1.z*coor1.z);
}

const double Link (Coordinate& k1, Coordinate& k2, Coordinate& r)
{
	double norm_1, norm_2, norm_temp, i, j, k;
	Coordinate norm_a, norm_b;

        norm_a = Cross(r, k1);
        norm_b = Cross(r, k2);

        norm_1 = Size(norm_a);
        norm_2 = Size(norm_b);

	if(norm_1 <= ZERO || norm_2 <= ZERO)  
	{
	    return 0;
        }
        else  
	{
	    norm_temp = norm_1 * norm_1 * norm_2 * norm_2;
            i = norm_b * norm_a / norm_temp;
            j = -1 * Size(r) * (r * Cross(k1, k2)) / norm_temp;
            k = atan2(j, -i);
            return k;
         }

}

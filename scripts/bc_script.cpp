#include <iostream>
#include <fstream>
#include <string>

using namespace std;


int main(int argc, char** argv)
{
    fstream newfile1;
    fstream newfile2;

    double count_file1 = 0;
    double count_file2 = 0;
    double count_both = 0;
    int file2_flag = 0;
    int end_of_file1 = 0;
    int end_of_file2 = 0;
    std::string kmer1;
    std::string kmer2;
    double count1;
    double count2;

    newfile1.open(argv[1],ios::in);
    newfile2.open(argv[2],ios::in);

    string tp;

    while(!end_of_file1)
    {
        if(getline(newfile1, tp))
        {
            size_t pos = tp.find('\t');
            kmer1 = tp.substr(0,pos);
            count1 = stoi(tp.substr(pos, tp.length() - pos));
            count_file1 += count1;
        }
        else
        {
            end_of_file1 = 1;
            break;
        }
        
        while(!end_of_file2)
        {
            if(file2_flag == 0)
            {
                if(getline(newfile2, tp))
                {
                    size_t pos = tp.find('\t');
                    kmer2 = tp.substr(0,pos);
                    count2 = stoi(tp.substr(pos, tp.length() - pos));   
                    count_file2 += count2;
                }      
                else
                {
                    end_of_file2 = 1;
                    break; 
                }                
            }
            
            if(kmer1.compare(kmer2) > 0)
            {
                file2_flag = 0;
                continue;

            }
            else if(kmer1.compare(kmer2) == 0)
            {
                count_both += min(count1, count2);
                file2_flag = 0;
                break;
            }
            else
            {
                file2_flag = 1;
                break;
            }

        }

        if(end_of_file2 == 1)
        {
            break;
        }

    }

    while(end_of_file1 == 0 && end_of_file2 == 1 && getline(newfile1, tp))
    {
        size_t pos = tp.find('\t');
        kmer1 = tp.substr(0,pos);
        count1 = stoi(tp.substr(pos, tp.length() - pos));
        count_file1 += count1;
    }

    while(end_of_file1 == 1 && end_of_file2 == 0 && getline(newfile2, tp))
    {
        size_t pos = tp.find('\t');
        kmer2 = tp.substr(0,pos);
        count2 = stoi(tp.substr(pos, tp.length() - pos));
        count_file2 += count2;
    }

    newfile1.close();
    newfile2.close();

    double dist = 1.0 - (2.0 * ((count_both)/(count_file1 + count_file2)));

    std::cout << "COMMON IN BOTH DATASETS: " << count_both << std::endl;
    std::cout << "TOTAL IN THE TWO DATASETS: " << count_file1 << " " << count_file2 << std::endl;
    std::cout << "BC DISSIMILARITY: " << dist << std::endl;

    return 0;

}

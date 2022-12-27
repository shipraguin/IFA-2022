#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <chrono>
using namespace std::chrono;

std::vector<int> suftab (std::string reference) {

    std::vector<std::string> suffixe (reference.size());
    std::vector<int> suffixe_int (reference.size()); 
    int i;
    std::string temp;

    for (i = 0; i < reference.size(); i++){
         suffixe[i] = reference.substr(i, reference.size()-i);
    }
    for (int i = 0; i < reference.size(); i++){
         for (int j = i + 1; j < reference.size(); j++){
            if (suffixe[i] > suffixe[j]){
                 temp = suffixe[i];
                 suffixe[i] = suffixe[j];
                 suffixe[j] = temp;
            }
         }
    }
    // Next we create the proper suffix array consisting of integers which
    // we return as the only output.
    for (int i = 0; i < reference.size(); i++){
        for (int j = 0; j < reference.size(); j++){
        if (suffixe[i] == reference.substr(j, reference.size() - j)){
            suffixe_int[i] = j;
        }
      }
    }
    return suffixe_int ;
}
int suffix_array_Lp(std::string p, std::vector<int> suftab, int array_length, 
        const std::string& reference) {

        int min_0 = std::min(p.size(), reference.size()-suftab[0] + 1);// +1 added later
        int min_l1 = std::min(p.size(), reference.size()- 
        suftab[array_length-1] + 1);
        int n = reference.size() - 2; // the index of last character

        if(p.compare(reference.substr(suftab[0], min_0))<= 0){
        return 0;
        } else if(p.compare(reference.substr(suftab[array_length-1],min_l1)) > 0){
        return n + 2; 
       }
      else {
        int L = 0;
        int R = array_length - 1; 
        while (R-L > 1){
            int M = ceil((L+R)/2);
            int min_M1 =  std::min(p.size(), reference.size() - suftab[M] + 1);
            if (p.compare(reference.substr(suftab[M], min_M1)) <= 0){ 
                R = M;
            }
            else{
                L = M;
           }
        } 
        return R;
    }
}
// We proceed similarly for the lower bound of the pattern we are looking for 
// to find in the reference by using binary search. Whenever the pattern is
// lesser than the suffixe at the midpoint of the suffixe array, we move upward
// the prefixed value of R, such that its values never goes beyond the pattern itself.
// We also take care of extreme cases that may appear and account for cases in which
// the reads may contain other nucleotides than the exected ones.

int suffix_array_Rp(int Lp, std::string p, std::vector<int> suftab, int array_length, const std::string& reference) {
    int min_l2 = std::min(p.size(), reference.size()-suftab[array_length-1] + 1);
    if (Lp == array_length - 1){
        return Lp;
    }

    // if p is situated on the last suffix of the sufftab, i.e. the sorted s.array
    if (p.compare(reference.substr(suftab[array_length-1],min_l2))==0 ){
         return array_length - 1; 
        }  
    else {
        int L = Lp;
        int R = array_length-1;
        while (R-L > 1){ 
            int M = ceil((L+R)/2);
            int min_M2 = std::min(p.size(), reference.size() - suftab[M] + 1);
            int min_M3 = std::min(p.size(), reference.size() - suftab[M+1] + 1);
            
            if (p.compare(reference.substr(suftab[M], min_M2)) < 0){
                R = M;
            }
            else if ((p.compare(reference.substr(suftab[M],
            min_M2))==0)
            && (p.compare(reference.substr(suftab[M+1],min_M3)) < 0)){
                R = M;
            }
            else {
                L = M;
            }          
        }
        if (p.compare(reference.substr(suftab[R], min_l2))==0){
            return R;
            }
        else {
            return R-1;
            }    
    }
}
// After having found the upper and lower bound in the suffix array where the pattern
// is found, we calculate the total number of the occurrence of the pattern by the 
// formula Rp-Lp+1.
int count(std::string p,std::vector<int> suftab, int array_length, const std::string& reference){
    int Lp = suffix_array_Lp(p, suftab, array_length, reference); 
    int Rp = suffix_array_Rp(Lp, p, suftab, array_length, reference);
    if(Lp == Rp && (p.compare(reference.substr(suftab[Rp],
            std::min(p.size(), reference.size()-suftab[array_length-1] + 1))) !=0)){
                return 0;
            }

    if (Lp > array_length - 1){
        return 0;
    }
    else if(p.compare(reference.substr(suftab[1], 
            std::min(p.size(), reference.size() - suftab[1] + 1))) < 0){
        return 0;
            } 
    else{    
        return Rp - Lp + 1;
    }
}
// The main function through which we provide the required output of the
// frequency of the pattern in the reference text.
int main() {
    // We read the reference file as a single line
    std::ifstream infile("hg38_10k.txt");
    std::string buffer; getline(infile, buffer); 
    int size = buffer.size();
    std::cout<<"The size of the big file: "<< size << std::endl;

    // We read the file containg the reads in a line-by-line manner
    //int array_counts[1000]; // change this to 1000, 10000, 97999
    //int size_array =  sizeof(array_counts) / sizeof(int);
    int sum = 0;
    std::ifstream input("illumina_reads_100.txt");

    // We call the suffix array generating function a single time and
    // measure the time needed to create it
    std::vector<int> suffix_array = suftab(buffer);
    auto start = high_resolution_clock::now();
    for (int i = 0; i < 1000; i++){
        std::string line; getline(input, line); 
        int n = count(line, suffix_array, size, buffer);
        sum = sum + n;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time execution for binary search: " << duration.count() << " microseconds." << std::endl;

    std::cout << "Total count: " << sum << "\n";
  return 0;
}

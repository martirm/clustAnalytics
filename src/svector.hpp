#ifndef SVECTOR_HPP // include guard
#define SVECTOR_HPP

#include <map>
#include <vector>
#include <random>
#include <Rcpp.h>


template <typename T>
class SVector{
    /* Allows sampling elements
     * Implemented as a unordered_map to keep the cost of insertions and deletions constant
     * std::sample function might have linear running time, which would be a problem (look into this)
     */
        std::vector<T> v;
        std::map<T, int> m;
        int l;
        std::mt19937 generator;
        std::uniform_int_distribution<int> distribution;
    public:
        SVector(){
            v = std::vector<T>();
            l = 0;
            m = std::map<T, int>();
            std::random_device device;
            generator = std::mt19937(device());
            distribution = std::uniform_int_distribution<int> (0, l-1);
        }
        SVector(std::vector<T> v0){
            v = v0;
            l = v.size();
            for (int i=0; i<l; ++i){
                m[v[i]] = i;
            }
            std::random_device device;
            generator = std::mt19937(device());
            distribution = std::uniform_int_distribution<int> (0, l-1);
        }
        void update_distribution(){
            distribution = std::uniform_int_distribution<int> (0, l-1);
        }
        void insert(T a){
            //does nothing if a is already there
            auto it = m.find(a);
            if (it == m.end()){
                v.push_back(a);
                m[a] = v.size()-1;
                ++l;
                update_distribution();
            }
        }
        void remove(T a){
            auto it = m.find(a);
            if (it != m.end()){
                int a_pos = it->second;
                m.erase(it);
                T last = v.back();
                if (a != last){
                    v[a_pos] = last;
                    m[last] = a_pos;
                }
                v.pop_back();
                --l;
                update_distribution();
            }
        }
        
        bool belongs(T a){
            // returns true if a belongs to the vector, false otherwise
            auto it = m.find(a);
            return not (it == m.end());
        }
        
        T rand_el_old(){
            int i = distribution(generator);
            std::cout << "random index: " << i << ", l=" << l << std::endl;
            return v[i];
            //return (v[distribution(generator)]);
        }
        T rand_el(){
            auto i = Rcpp::sample(l, 1);
            return v[i[0]-1];
        }
        // vector<T> sample(int n){
        //     std::vector<T> s(n);
        //     std::sample(v.begin(), v.end(), s.begin(), n);
        //     return (s);
        // }
};


#endif /* SVECTOR_HPP */

#ifndef CANTOR_HASH_H // include guard
#define CANTOR_HASH_H

class CantorHash {
    public:
        std::size_t operator()(const std::pair<int,int>& p) const {
            return (p.first + p.second) * (p.first + p.second + 1) / 2 + p.second;
        }
};


#endif /*CANTOR_HASH_H*/
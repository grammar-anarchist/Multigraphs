#ifndef MULTIGRAPHS_ZP_H
#define MULTIGRAPHS_ZP_H

class Zp {
public:
    int prime;
    int num;
    Zp(): prime(2), num(0) {}
    explicit Zp(int prime, int num): prime(prime), num(num) {}
    Zp &operator=(int a) {
        if (a >= 0) {
            num = a % prime;
        } else if (a % prime == 0) {
            num = 0;
        } else {
            num = a % prime + prime;
        }
        return *this;
    }

    // must have the same prime
    Zp operator*(Zp a) const {
        Zp ans(prime, this->num * a.num);
        return ans;
    }
    Zp operator+(Zp a) const {
        Zp ans(prime, this->num + a.num);
        return ans;
    }
    bool operator==(int a) const {
        return this->num == a % prime;
    }
    // a must have the same prime
    Zp operator/(Zp a) {
        int rev = 1;
        for (int i = 0; i != prime - 2; ++i) {
            rev *= a.num;
            rev %= prime;
        }
        auto other = Zp(prime, rev);
        return (*this) * other;
    }
};

#endif //MULTIGRAPHS_ZP_H

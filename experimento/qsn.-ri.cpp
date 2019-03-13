//-----------------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
//-----------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
typedef long double real;
//------------------------------------------------------------------------------
class clEvent{
   public:
      int  nq;
      real iat, st, at, bs, es, idt;
};
//------------------------------------------------------------------------------
class clQS{
   public:
      real l, m, T, U, E[7], V[7];
      vector<clEvent> Event;
};
//------------------------------------------------------------------------------
class clQSN{
   private:
      size_t s;
      vector<size_t> S;
      real F     ( real );
      void Open  ( size_t, size_t );
      void Close ( size_t, size_t );
      void Nq    ( void );

   public:
      vector<clQS> QS;

      void Config   ( int, vector<real>, int, vector<size_t> );
      void Simulate ( int );
      void Statistic( int );
};
//------------------------------------------------------------------------------
void clQSN::Config( int qs, vector<real> D, int s, vector<size_t> S ){
    this->s = s;
    this->S = S;
    QS.clear();
    for( int i = 0; i < qs; i++ ){
         clQS x;
         x.l = 0.0;
         if( i == 0 ) x.l = D[i];
         x.m = D[i+1];
         QS.push_back(x);
    }
}
real clQSN::F( real p ){
   real u = (rand()+1.0)/(RAND_MAX+2.0); // u in (0,1)
   return -p*log(u);
}
void clQSN::Open( size_t i, size_t f ){
     clEvent X, Xa;
     Xa = QS[i].Event[ QS[i].Event.size()-1 ];
     if( i == 0 ){
         X.iat = F( QS[i].l );
         X.st  = F( QS[i].m );
         X.at  = Xa.at + X.iat;
         X.bs  = X.at > Xa.es ? X.at : Xa.es;
         X.es  = X.bs + X.st;
         X.idt = X.es - Xa.es;
         QS[i].Event.push_back(X);
//         Nq(i);
     }
     Xa    = QS[i].Event[ QS[i].Event.size()-1 ];
     X.iat = Xa.idt;
     X.st  = F( QS[f].m );

     Xa    = QS[f].Event[ QS[f].Event.size()-1 ];
     X.at  = Xa.at + X.iat;
     X.bs  = X.at > Xa.es ? X.at : Xa.es;
     X.es  = X.bs + X.st;
     X.idt = X.es - Xa.es;
     QS[f].Event.push_back(X);
//     Nq(f);
}
void clQSN::Close( size_t i, size_t f ){
     clEvent X, Xa;

     Xa    = QS[i].Event[ QS[i].Event.size()-1 ];
     X.iat = Xa.idt;
     X.st  = F( QS[f].m );

     Xa = QS[f].Event[ QS[f].Event.size()-1 ];
     X.at  = Xa.at + X.iat;
     X.bs  = X.at > Xa.es ? X.at : Xa.es;
     X.es  = X.bs + X.st;
     X.idt = X.es - Xa.es;
     QS[f].Event.push_back(X);
//     Nq(f);
}
void clQSN::Nq( void ){
   for( size_t qs = 0; qs < QS.size(); qs++ ){
        for( size_t e = 1; e < QS[qs].Event.size(); e++ ){
             size_t c = e-1;
             QS[qs].Event[e].nq = 0;
             while( QS[qs].Event[e].at < QS[qs].Event[c].es ){
                    QS[qs].Event[e].nq += 1;
                    c--;
             }
        }
   }
}
void clQSN::Simulate( int N ){
   for( size_t i = 0; i < s; i++ )
        QS[ S[i] ].Event.push_back({0,0,0,0,0,0,0});
   for( int e = 1; e < N; e++ ){
        for( size_t i = 1; i < s; i++ ){
             Open(S[i-1],S[i]);
            //Close(S[i-1],S[i]);
        }
   }
}
void clQSN::Statistic( int Ni ){
   Nq();
   for( size_t qs = 0; qs < QS.size(); qs++ ){
        size_t N = QS[qs].Event.size();
        real   x[7], Sx[7], Sxx[7];
        for( size_t s = 0; s < 7; s++ )
             Sx[s] = Sxx[s] = 0.0;
        for( size_t e = Ni; e < N; e++ ){
             clEvent X = QS[qs].Event[e];
             x[0] = X.iat    ;
             x[1] = X.st     ;
             x[2] = X.nq     ;
             x[3] = X.idt    ;
             x[4] = X.bs-X.at; //w
             x[5] = X.es-X.bs; //s
             x[6] = X.es-X.at; //r
             for( int s = 0; s < 7; s++ ){
                  Sx [s] += x[s];
                  Sxx[s] += x[s]*x[s];
             }
        }
        QS[qs].T = QS[qs].Event[N-1].es-QS[qs].Event[Ni-1].bs;
        QS[qs].U = Sx[1]/QS[qs].T;
        for( int s = 0; s < 7; s++ ){
             QS[qs].E[s] = Sx [s]/(N-Ni);
             QS[qs].V[s] = Sxx[s]/(N-Ni)-QS[qs].E[s]*QS[qs].E[s];
        }
   }
}
string RI( void ){
   int R  = 10,
       N  = 1000,
       Ni = 0.9*N,
       qs = 5,
       s  = 17;
   vector<size_t> S = { 0, 1, 2, 1, 3, 1, 2, 1, 3, 1, 4, 1, 2, 1, 4, 1, 0 };
   vector<real>   D = { 0.0288, 0.0027, 0.002, 0.005, 0.0222, 0.0284 };

   clQSN QSN[R];
   stringstream str;

   str << fixed;
   str.precision(3);

   for( int r = 0; r < R; r++ ){
        cout << " Calculando RI: " << r << endl;
        srand( time(NULL)/(r+1) );
        QSN[r].Config(qs,D,s,S);
        QSN[r].Simulate(N);
        QSN[r].Statistic(Ni);
   }

   str << "<html><body>"
       << "<br>QSN Table"
       << "<br>QS Number:" << qs
       << "<br>Sequence:";

    for( int i = 1; i < s; i++ ){
         str << S[i-1] << "-" << S[i] << ";";
    }
    str << "</table>"
           "<br><br><table border='1' cellpadding='0' cellspacing='0'>"
           "<tr><td>QSN<td>l<td>m";
    for( int q = 0; q < qs; q++ ){
         str << "<tr><td>" << q;
         if( q == 0) str << "<td>" << D[q];
         else        str << "<td>-";
         str << "<td>" << D[q+1];
    }
    str << "</table>"
        << "<br><br><table border='1' cellpadding='0' cellspacing='0'>"
        << "<tr><td>r";
    for( int q = 0; q < qs; q++ ){
         str << "<td>QS"
             << "<td>T"
//             << "<td>E[iat]"
//             << "<td>E[st]"
             << "<td>E[nq]"
//             << "<td>E[idt]"
             << "<td>E[w]"
//             << "<td>E[s]"
//             << "<td>E[r]"
//           << "<td>V[iat]<td>V[st]<td>V[nq]<td>V[idt]<td>V[w]<td>V[s]<td>V[r]"
             << "<td>U";
//             << "<td>p0";
    }

    for( int r = 0; r < R; r++ ){
         str << "<tr><td>" << r+1;
         for( size_t qs = 0; qs < QSN[r].QS.size(); qs++ ){
              str << "<td>" << qs
                  << "<td>" << QSN[r].QS[qs].T
//                  << "<td>" << QSN[r].QS[qs].E[0]           // iat
//                  << "<td>" << QSN[r].QS[qs].E[1]           // st
                  << "<td>" << QSN[r].QS[qs].E[2]             // nq
//                  << "<td>" << QSN[r].QS[qs].E[3]           // idt
                  << "<td>" << QSN[r].QS[qs].E[4]             // w
//                  << "<td>" << QSN[r].QS[qs].E[5]           // s
//                  << "<td>" << QSN[r].QS[qs].E[6]           // r
//                  << "<td>" << QSN[r].QS[qs].V[0]
//                  << "<td>" << QSN[r].QS[qs].V[1]
//                  << "<td>" << QSN[r].QS[qs].V[2]
//                  << "<td>" << QSN[r].QS[qs].V[3]
//                  << "<td>" << QSN[r].QS[qs].V[4]
//                  << "<td>" << QSN[r].QS[qs].V[5]
//                  << "<td>" << QSN[r].QS[qs].V[6]
                  << "<td>" << QSN[r].QS[qs].U;
//                  << "<td>" << 1.0-QSN[r].QS[qs].U;
         }
    }
    return str.str();
}
int main( void ){
    string   str = RI();

    ofstream fo;
    fo.open("out.QSN.MM1-RI-C++.html");
    replace(str.begin(), str.end(),'.',',');
    fo << str;
    fo.close();

    return 0;
}
//------------------------------------------------------------------------------
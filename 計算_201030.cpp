#include <bits/stdc++.h>

#ifdef __FreeBSD__
#include <floatingpoint.h>               
#endif
#ifdef __FreeBSD__
   fpsetprec(FP_PE);                
#endif


#define rep(i,n) for(int i = 0; i < (n); ++i)
#define srep(i,s,t) for(int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;

const int m = 4; // 塩基の種類数(A,C,G,T)

/*
    パラメータR,qはmodel_params.txtより受け取る。(これは木立先生が最尤法により推定したもの)
    Rから固有値、固有ベクトルを計算しU,U_inv,lambdaを得る。(これはPythonで計算)
*/

// initial distribution q : 初期確率 (和は1)
vector<double> q = {0.2308288445620057, 0.2691711554379943, 0.2691711554379943, 0.2308288445620057}; 

// transition matrix P : 遷移確率行列 (P[i][j] : jからiに遷移する確率、よってj固定でiを動かした(各列の)和は1)
// 今回は連続時間モデルのため使用しない。
/*
double P[m][m]; 
double Pd[m][m]; // diagonal part of matrix P : Pの対角成分 (つまりP[i][i]:遷移(塩基置換)しない確率。またP[i][j](i!=j)の値は0)
double Ps[m][m]; // non diagonal part of matrix P : Pの非対角成分 (つまり遷移(塩基置換)する確率。またP[i][i]の値は0)
*/

vector<vector<double>> R = {
    {-1.0439782893776948,	0.1585886443107131,	0.5657877891067001,	0.19927682168764793},
    {0.25488550029885526,	-0.96229058720693,	0.23791415378951666,	0.5898159673911916},
    {0.5898159673911916,	0.23791415378951666,	-0.9622905872069297,	0.25488550029885526},
    {0.19927682168764793,	0.5657877891067001,	0.1585886443107131,	-1.0439782893776948}
};

vector<vector<double>> U = {
    {0.460305948054125, 0.509758760460158, -0.500000000000005, -0.537282769914936}, 
    {0.536766647795849, -0.490046942786223, 0.499999999999995, -0.459703409985757}, 
    {0.536766647795849, 0.490046942786223, 0.500000000000005, 0.459703409985748}, 
    {0.460305948054125, -0.509758760460158, -0.499999999999995, 0.537282769914945} 
};

vector<vector<double>> U_inv = {
    {0.501467999502849, 0.501467999502850, 0.501467999502850, 0.501467999502850}, 
    {0.461891278801117, -0.539839862579025, 0.539839862579025, -0.461891278801116}, 
    {-0.538342594140075, 0.461657405859925, 0.461657405859934, -0.538342594140066}, 
    {-0.492379225734085, -0.512184857963776, 0.512184857963767, 0.492379225734095} 
};

vector<vector<double>> lambda = {
    {0, 0, 0, 0},
    {0, -0.851801914984303, 0, 0},
    {0, 0, -1.569077901107460, 0},
    {0, 0, 0, -1.591657937077486}
};

vector<vector<double>> Rd(m, vector<double>(m));
vector<vector<double>> Rs(m, vector<double>(m));

/* ここまでは与えられる。 */

vector<vector<double>> dot(vector<vector<double>> a, vector<vector<double>> b) {
    vector<vector<double>> res(m, vector<double>(m));
    rep(i,m) rep(j,m) res[i][j] = 0;
    rep(i,m) {
        rep(j,m) {
            rep(k,m) {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return res;
}

vector<vector<double>> calcPt(double t) {
    vector<vector<double>> et(m, vector<double>(m));
    rep(i,m) {
        rep(j,m) {
            if(i == j) {
                et[i][j] = exp(t*lambda[i][i]);
            }else {
                et[i][j] = 0;
            }
        }
    }
    return dot(dot(U, et), U_inv);
}

double calcK(double a, double b) {
    if(a != b) {
        return (exp(a) - exp(b)) / (a - b);
    }else {
        return exp(a);
    }
}


struct node {
    int id; // 頂点番号(今のところ使っていない)
    int parent; // 親の頂点番号(rootの場合は-1)
    int sonLeft; // 左側の子の頂点番号(leafの場合は-1)
    int sonRight; // 右側の子の頂点番号(leafの場合は-1)
    int flagLeaf; // 葉か否か
    int flagRoot; // 根か否か
    int flagHumanPath; // (頂点が)ヒトのパス上に属しているかどうか
    double distParent; // 親までの枝の距離
    double distLeft; // 左側の子までの距離
    double distRight; // 右側の子までの距離
    int base; // 塩基(leafのみ、番号)
    string name; // 生物種名(leaf以外は空)
    int lca; // ヒトとの最近共通祖先の頂点番号
};


// パーサー
vector<node> parseNewick(string s) {
    int n = s.size();

    stack<int> st, st2;
    int cnt = 0;

    int now = 0;
    int lr = 0;

    vector<node> v;

    rep(i,n) {
        // cout << i << ' ' << cnt << ' ' << ' ' << s[i] << ' ' << now << endl;
        if(s[i] == ';') {
            break;
        }else if(s[i] == '(') {
            node nd;
            if(st.size() == 0) {
                nd.parent = -1;
            }else {
                nd.parent = st.top();
            }
            st.push(cnt);
            st2.push(lr);
            now = cnt;
            cnt++;
            nd.flagLeaf = 0;
            v.push_back(nd);
            lr = 0;
        }else if(s[i] == ')') {
            st.pop();
            lr = st2.top();
            st2.pop();
        }else if(s[i] == ':') {
            string num;
            while(s[i+1] == '.' || ('0'<=s[i+1]&&s[i+1]<='9')) {
                num += s[i+1];
                i++;
            }
            now = st.top();
            // cout << "num " << num << ' ' << now << endl;
            if(lr == 0) {
                v[now].distLeft = stod(num);
            }else {
                v[now].distRight = stod(num);
            }
        }else if(s[i] == ',') {
            lr = 1;
        }else {
            string name;
            name += s[i];
            while(s[i+1] != ':') {
                name += s[i+1];
                i++;
            }
            node nd;
            nd.flagLeaf = 1;
            nd.sonLeft = -1;
            nd.sonRight = -1;
            nd.distLeft = 0;
            nd.distRight = 0;
            nd.name = name;
            nd.parent = now;
            v.push_back(nd);
            cnt++;
        }
    }

    // flagRoot埋め
    rep(i,cnt) {
        if(v[i].parent == -1) v[i].flagRoot = 1;
        else v[i].flagRoot = 0;
    }

    // sonLeft, sonRight, distParent埋め
    /*
        頂点番号が小さいほうが必ずLeft側なので無理矢理そうする。
    */
    vector<int> sons[cnt];
    rep(i,cnt) {
        if(v[i].flagRoot) continue;
        sons[v[i].parent].push_back(i);
    }
    rep(i,cnt) {
        sort(sons[i].begin(), sons[i].end());
        if(sons[i].size() == 2) {
            v[i].sonLeft = sons[i][0];
            v[sons[i][0]].distParent = v[i].distLeft;
            v[i].sonRight = sons[i][1];
            v[sons[i][1]].distParent = v[i].distRight;
        }
    } 

    

    
    /*
        ここで、頂点番号が小さい方を左側の子とする。
        今回は、これでmafファイルと整合性が取れそう。
        ！！！バグ！！！
        上のパースはDFSでやってるけどこの変数埋めはBFSでやってて頂点番号がずれてる。
        なので、削除
    */
    /*
    queue<int> que;
    rep(i,cnt) {
        if(v[i].flagLeaf) que.push(i);
    }
    
    int sonCount[cnt] = {};

    while(!que.empty()) {
        int x = que.front();
        que.pop();
        if(v[x].parent == -1) continue;
        int px = v[x].parent;
        if(sonCount[px] == 0) {
            v[px].sonLeft = x;
            v[x].distParent = v[px].distLeft;
            sonCount[px]++;
        }else {
            v[px].sonRight = x;
            v[x].distParent = v[px].distRight;
            que.push(px);
        }
    }
    */


    // flagHumanPath埋め
    rep(i,cnt) v[i].flagHumanPath = 0;
    rep(i,cnt) {
        if(v[i].name == "hg38") {
            v[i].flagHumanPath = 1;
            int now = i;
            while(v[now].parent != -1) {
                now = v[now].parent;
                v[now].flagHumanPath = 1;
            }
            break;
        }
    }


    // 最近共通祖先埋め
    rep(i,cnt) {
        int now = i;
        while(!v[now].flagHumanPath) now = v[now].parent;
        v[i].lca = now;
    }

    return v;
}

vector<vector<double>> calcQ(node v) {
    vector<vector<double>> res(m, vector<double>(m));
    rep(i,m) rep(j,m) res[i][j] = 0;

    double t = v.distParent;

    rep(i,m) {
        rep(j,m) {
            rep(l,m) {
                rep(h,m) {
                    res[i][j] += t * Rs[i][l] * U[l][h] * U_inv[h][j] * calcK(t*Rd[i][i], t*lambda[h][h]);
                }
            }
        }
    }

    return res;
}

string allTissueName[41] = {
    "adipose_tissue",
    "blood",
    "blood_vessel",
    "brain",
    "esophagus",
    "eye",
    "female_gonad",
    "gallbladder",
    "heart",
    "internal_male_genitalia",
    "kidney",
    "large_intestine",
    "liver",
    "lung",
    "lymph_node",
    "meninx",
    "olfactory_region",
    "pancreas",
    "parotid_gland",
    "penis",
    "placenta",
    "prostate_gland",
    "salivary_gland",
    "skeletal_muscle_tissue",
    "skin_of_body",
    "small_intestine"      ,
    "smooth_muscle_tissue"   ,
    "spinal_cord",
    "spleen",
    "stomach",
    "submandibular_gland",
    "testis",
    "throat",
    "thymus",
    "thyroid_gland",
    "tongue",
    "tonsil",
    "umbilical_cord",
    "urinary_bladder",
    "uterus",
    "vagina"
};


void RdRsInIt() {
    // Rd, Rs初期化 
    rep(i,m) {
        rep(j,m) {
            if(i == j) {
                Rd[i][j] = R[i][j];
                Rs[i][j] = 0;
            }else {
                Rd[i][j] = 0;
                Rs[i][j] = R[i][j];
            }
        }
    }
}

vector<int> leafNodeNumbersInIt(const vector<node>& V) {
    int n = V.size();
    vector<int> leafNodeNumbers;
    rep(i,n) {
        if(V[i].flagLeaf) {
            leafNodeNumbers.push_back(i);   
        }
    }
    return leafNodeNumbers;
}

vector<vector<double>> nuInIt(const vector<node>& V) {
    int n = V.size();
    vector<vector<double>> nu(n, vector<double>(m));
    rep(i,n) {
        if(V[i].flagLeaf) {
            if(V[i].base == -1) {
                rep(j,m) {
                    nu[i][j] = 1;
                }
            }else {
                rep(j,m) {
                    if(j == V[i].base) nu[i][j] = 1;
                    else nu[i][j] = 0;
                }
            }
        }else {
            rep(j,m) nu[i][j] = 0;
        }
    }
    return nu;
}


vector<vector<double>> alphaCalc(const vector<node>& V, const vector<vector<double>>& nu, const vector<vector<vector<double>>>& Pt) {
    int n = V.size();
    vector<vector<double>> alpha(n, vector<double>(m));

    queue<int> que;
    rep(i,n) {
        if(V[i].flagLeaf) {
            que.push(i);
        }
    }
    int cntAlphaSon[n] = {};
    while(!que.empty()) {
        int x = que.front();
        que.pop();
        // alpha[x]計算
        if(V[x].flagLeaf) {
            rep(j,m) alpha[x][j] = nu[x][j];
        }else {
            int left = V[x].sonLeft;
            int right = V[x].sonRight;
            double tmpLeft[m] = {};
            double tmpRight[m] = {};
            rep(i,m) {
                rep(j,m) {
                    tmpLeft[i] += alpha[left][j] * Pt[left][j][i];
                    tmpRight[i] += alpha[right][j] * Pt[right][j][i];
                }
            }
            rep(i,m) alpha[x][i] = tmpLeft[i] * tmpRight[i];
        }

        // 親に伝播
        if(V[x].flagRoot) continue;
        cntAlphaSon[V[x].parent]++;
        if(cntAlphaSon[V[x].parent] == 2) que.push(V[x].parent);
    }

    return alpha;
}


vector<vector<double>> betaCalc(const vector<node>& V, const vector<vector<vector<double>>>& Pt, const vector<vector<double>>& alpha) {
    int n = V.size();
    vector<vector<double>> beta(n, vector<double>(m));

    int rootNumber = 0;
    rep(i,n) {
        if(V[i].flagRoot) rootNumber = i;
    }
    
    queue<int> que;
    que.push(rootNumber); // queは再利用

    while(!que.empty()) {
        int x = que.front();
        que.pop();
        if(V[x].flagLeaf) continue;
        // beta[x]計算
        if(V[x].flagRoot) {
            double tmpSon[m] = {};
            int son = 0;
            if(V[V[x].sonLeft].flagHumanPath) son = V[x].sonRight;
            else son = V[x].sonLeft;
            rep(i,m) {
                rep(j,m) {
                    tmpSon[i] += alpha[son][j] * Pt[son][j][i];
                }
            }
            rep(i,m) beta[x][i] = tmpSon[i] * q[i];
        }else {
            double tmpSon[m] = {};
            int son = 0;
            if(V[V[x].sonLeft].flagHumanPath) son = V[x].sonRight;
            else son = V[x].sonLeft;
            rep(i,m) {
                rep(j,m) {
                    tmpSon[i] += alpha[son][j] * Pt[son][j][i];
                }
            }

            double tmpParent[m] = {};
            int par = V[x].parent;
            rep(i,m) {
                rep(j,m) {
                    tmpParent[i] += Pt[x][i][j] * beta[par][j];
                }
            }
            rep(i,m) beta[x][i] = tmpSon[i] * tmpParent[i];
        }

        // ヒトのパスのみ辿る
        if(V[x].flagLeaf) continue;
        if(V[V[x].sonLeft].flagHumanPath) que.push(V[x].sonLeft);
        if(V[V[x].sonRight].flagHumanPath) que.push(V[x].sonRight);
    }

    return beta;
}


vector<vector<double>> gammaCalc(const vector<node>& V, const vector<vector<double>>& nu, const vector<vector<vector<double>>>& Pt, const vector<vector<double>>& alpha) {
    int n = V.size();
    vector<vector<double>> gamma(n, vector<double>(m));

    int humanNumber = 18; // ヒトの頂点番号は18
    rep(i,n) if(V[i].flagHumanPath && V[i].flagLeaf) humanNumber = i;

    queue<int> que;
    que.push(humanNumber); // queを再利用

    while(!que.empty()) {
        int x = que.front();
        que.pop();
        // gamma[x]計算
        if(V[x].flagLeaf) {
            rep(i,m) gamma[x][i] = nu[x][i];
        }else {
            int human, nonHuman;
            if(V[V[x].sonLeft].flagHumanPath) {
                human = V[x].sonLeft;
                nonHuman = V[x].sonRight;
            }else {
                human = V[x].sonRight;
                nonHuman = V[x].sonLeft;
            }
            double tmpHuman[m] = {};
            double tmpNonHuman[m] = {};
            double t = V[human].distParent;
            rep(i,m) {
                rep(j,m) {
                    if(j == i) tmpHuman[i] += gamma[human][j] * exp(t*Rd[j][i]);
                    tmpNonHuman[i] += alpha[nonHuman][j] * Pt[nonHuman][j][i];
                }
            }
            rep(i,m) gamma[x][i] = tmpHuman[i] * tmpNonHuman[i];
        }

        // 親を辿る
        if(V[x].flagRoot) continue;
        que.push(V[x].parent);
    }


    return gamma;
}


double calcZ(const vector<node>& V, const vector<vector<double>>& alpha) {
    double Z = 0;
    int n = V.size();
    int rootNumber = 0;
    rep(i,n) {
        if(V[i].flagRoot) rootNumber = i;
    }
    Z = 0;
    rep(i,m) Z += alpha[rootNumber][i] * q[i];
    return Z;
}


vector<double> calcProb(const vector<node>& V, const vector<vector<vector<double>>>& Q, const vector<vector<double>>& beta, const vector<vector<double>>& gamma, const double Z) {
    int n = V.size();
    vector<double> prob(n);
    rep(i,n) {
        prob[i] = 0;
        if(!V[i].flagHumanPath) continue;
        if(V[i].flagRoot) {
            rep(j,m) {
                prob[i] += gamma[i][j] * q[j];
            }
        }else {
            rep(j,m) {
                rep(k,m) {
                    prob[i] += gamma[i][j] * Q[i][j][k] * beta[V[i].parent][k];
                }
            }
        }
        prob[i] /= Z;
    }

    return prob;
}

void writeToOfstream(const string enhancerName, const int tissueCount, const vector<string>& tissues, const vector<vector<double>>& allProbabilities, ofstream& ofs) {
    ofs << enhancerName << endl;
    ofs << tissueCount << endl;
    rep(i,tissueCount) ofs << tissues[i] << endl;

    ofs << allProbabilities.size() << endl;
    rep(i,allProbabilities.size()) {
        rep(j,allProbabilities[i].size()) {
            if(j == allProbabilities[i].size() - 1) ofs << allProbabilities[i][j] << endl;
            else ofs << allProbabilities[i][j] << '\t';
        }
    }
}

void writeLengthOfHumanPathEdges(const vector<node>& V) {
    int n = V.size();
    vector<double> lengthOfEdgeToParent(19, 0);
    rep(i,n){
        if(V[i].flagHumanPath && i > 0){
            lengthOfEdgeToParent[18-i] = V[i].distParent;
        }
    }
    drep(i,19){
        // lengthOfEdgeToParent[i] /= lengthOfEdgeToParent[0]; // Homoの枝の長さで割る(Homoの枝の長さを1とする)、またrootの枝の長さは0になる
    } 
    rep(i,19) cout << lengthOfEdgeToParent[i] << endl;

    string fileNameOfs = "./lengthOfHumanPathEdges.txt";
    const char* cstrOfs = fileNameOfs.c_str();
    ofstream ofs(cstrOfs);
    rep(i,19) ofs << lengthOfEdgeToParent[i] << '\t';
    ofs << endl;
    ofs.close();
}

int group[100] = {
    18,
    17,
    16,
    15,
    14,
    13,13,13,13,
    12,12,
    11,
    10,
    9,9,9,9,9,9,9,9,9,9,9,9,9,
    8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
    7,7,7,7,7,7,7,
    6,6,6,
    5,
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
    3,
    2,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    0
};

// array jobで連番を渡すためコマンドライン引数を用いる。
int main(int argc, char* argv[]) {
    clock_t start,end;
	start = clock();

    RdRsInIt();

    /*
        頂点数および各頂点の情報はnewick formatからパースする。
        各葉の塩基はmafファイルより得る。
        各頂点の生物種名はconcestors.txtを参照する。
        
    */

    // newick format とりあえずハードコード
    string strOfNewick = "((((((((((((((((((hg38:0.00775018,panTro4:0.00822868):0.0026527499999999997,gorGor3:0.010601039999999999):0.011252,ponAbe2:0.02215023):0.0043521,nomLeu3:0.02739288):0.01356963,(((rheMac3:0.00436145,macFas5:0.00338455):0.00626117,papAnu2:0.00969236):0.005104030000000001,chlSab2:0.01492744):0.03065243):0.02518389,(calJac3:0.04132332,saiBol1:0.0397847):0.044187019999999994):0.07755081999999999,otoGar3:0.18058727000000002):0.02029684,tupChi1:0.21400775):0.00644765,(((speTri2:0.1690805,(jacJac1:0.20809766,((micOch1:0.11519311,(criGri1:0.05933852,mesAur1:0.07380833):0.03884194):0.030547840000000003,(mm10:0.08832276,rn6:0.09438034000000001):0.06934151):0.13350729):0.0666497):0.01082311,(hetGla2:0.10586679,(cavPor3:0.13521471999999998,(chiLan1:0.09233307,octDeg1:0.13478029):0.01750058):0.030002749999999998):0.12135533999999999):0.034000129999999996,(oryCun2:0.12924065,ochPri3:0.21915592):0.11766914):0.01655578):0.0269221,(((susScr3:0.14374649,((vicPac2:0.0208641,camFer1:0.018902509999999997):0.11137344,((turTru2:0.00865893,orcOrc1:0.0070627):0.06986307,(panHod1:0.02135274,(bosTau8:0.05796649,(oviAri3:0.01748781,capHir1:0.01591031):0.008752919999999999):0.00193437):0.1271921):0.02537756):0.00504363):0.050787259999999994,(((equCab2:0.09226394,cerSim1:0.07078421):0.04012774,(felCat8:0.10474937000000001,(canFam3:0.10909532000000001,(musFur1:0.10936401,(ailMel1:0.07146708,(odoRosDiv1:0.03120006,lepWed1:0.02829569):0.03272742):0.0057826):0.02359607):0.02398091):0.06117955):0.005737900000000001,((pteAle1:0.00947624,pteVam1:0.011200410000000001):0.12212452000000001,(eptFus1:0.047568660000000006,(myoDav1:0.0338942,myoLuc2:0.02075108):0.02777062):0.11727507):0.03916131):0.00477061):0.01457266,(eriEur2:0.27863752999999997,(sorAra2:0.32450049999999997,conCri1:0.1888113):0.01975367):0.04373776):0.02739492):0.02383223,(((((loxAfr3:0.09977239,eleEdw1:0.25524171):0.00410012,triMan1:0.08699928):0.020550679999999998,(chrAsi1:0.16558979000000001,echTel2:0.25585381):0.0265092):0.00467971,oryAfe1:0.13653174):0.06693795,dasNov3:0.18667159000000003):0.012414650000000001):0.25944848,(monDom5:0.15918422,(sarHar1:0.15366878,macEug2:0.16562339999999998):0.03215452):0.27158742):0.06874389,ornAna1:0.55265026):0.14105729,(((((colLiv1:0.16760753,((falChe1:0.017144,falPer1:0.0159901):0.11970699,(((ficAlb2:0.09557994,((zonAlb1:0.07710477,geoFor1:0.048785800000000004):0.04279451,taeGut2:0.08683674):0.024226679999999997):0.015395399999999998,pseHum1:0.08371719999999999):0.15001801,(melUnd1:0.07659597,(amaVit1:0.06616203999999999,araMac1:0.08300262):0.04501083):0.08582126000000001):0.016193310000000002):0.00859234):0.03936694,(anaPla1:0.13819359,(galGal4:0.06739205,melGal1:0.0839683):0.11293024000000002):0.04598216):0.23027509999999998,allMis1:0.25588376):0.032607819999999996,((cheMyd1:0.06652824,chrPic2:0.06479837):0.047037899999999994,(pelSin1:0.05088364,apaSpi1:0.10378230000000001):0.11141492):0.11334820000000001):0.08739725,anoCar2:0.55254463):0.13859923999999998):0.10800258,xenTro7:0.90575609):0.04581508,latCha1:0.63016434):0.19984434,(((((((tetNig2:0.23367285,(fr3:0.032460420000000004,takFla1:0.07895543):0.17435836999999998):0.27523667,(oreNil2:0.049694,(neoBri1:0.05999335,(hapBur1:0.0324577,(mayZeb1:0.02476171,punNye1:0.02793736):0.00861032):0.019635720000000002):0.026521680000000002):0.27501105):0.00833543,(oryLat2:0.3946054,xipMac1:0.36460464000000004):0.09294102):0.04220054,gasAcu1:0.32457885):0.17418677999999999,gadMor1:0.58870043):0.16872787,(danRer10:0.54169416,astMex1:0.45580021000000004):0.23325083):0.17763958999999999,lepOcu1:0.56466518):0.16834161000000003):0.09985077,petMar2:1.14343783);";
    vector<node> V = parseNewick(strOfNewick); // 頂点配列
    int n = V.size(); // 頂点数
    vector<int> leafNodeNumbers = leafNodeNumbersInIt(V); // 生物種の頂点番号
    int numOfLeaf = leafNodeNumbers.size(); // 生物種数(100)  
    cout << "To parse newick and convert verticle vector finished." << endl;

    // writeLengthOfHumanPathEdges(V);
    
    // mafファイル(をパースしたもの)読み込み
    string fileNameIfs = "./enhancers_out/enhancers_out_";
    string tmpArgv = argv[1];
    if(tmpArgv.size() == 1) tmpArgv = '0' + tmpArgv;
    fileNameIfs = fileNameIfs + tmpArgv + ".txt";
    const char* cstrIfs = fileNameIfs.c_str();
    ifstream ifs(cstrIfs);

    string fileNameOfs = "./probabilities_out/probabilities_out_";
    tmpArgv = argv[1];
    if(tmpArgv.size() == 1) tmpArgv = '0' + tmpArgv;
    fileNameOfs = fileNameOfs + tmpArgv + ".txt";
    const char* cstrOfs = fileNameOfs.c_str();
    ofstream ofs(cstrOfs);

    
    
    int enhancerIte = 0;
    string line;

    int countTestis = 0, countEye = 0, countBrain = 0;
    double averageTestis = 0, averageEye = 0, averageBrain = 0;

    while(getline(ifs, line)) {
        string enhancerName = line;

        getline(ifs, line);
        int tissueCount = stoi(line);
        vector<string> tissues;
        int flagTestis = 1;
        int flagEye = 0;
        int flagBrain = 0;
        rep(i,tissueCount) {
            string temp;
            getline(ifs, temp);
            tissues.push_back(temp);
            if(temp == "testis") flagTestis = 1;
            if(temp == "eye") flagEye = 1;
            if(temp == "brain") flagBrain = 1;
        }

        vector<string> enhanserSequences(numOfLeaf);
        rep(i,numOfLeaf) {
            getline(ifs, enhanserSequences[i]);
        }

        vector<int> argmaxProb;
        vector<double> maxProbScore;
        vector<vector<double>> allProbabilities;

        rep(loop,enhanserSequences[0].size()) {
            string cKeep;
            int flagHumanBase = 0;
            int maxUnGap = 0;
            rep(i,numOfLeaf) {
                char c = enhanserSequences[i][loop];
                c = toupper(c);
                if(c == 'A') V[leafNodeNumbers[i]].base = 0;
                else if(c == 'C') V[leafNodeNumbers[i]].base = 1;
                else if(c == 'G') V[leafNodeNumbers[i]].base = 2;
                else if(c == 'T') V[leafNodeNumbers[i]].base = 3;
                else  V[leafNodeNumbers[i]].base = -1;
                cKeep += c;
                if(c=='A'||c=='C'||c=='G'||c=='T') maxUnGap = i;
                if(V[leafNodeNumbers[i]].flagHumanPath) {
                    if(V[leafNodeNumbers[i]].base == -1) break;
                    else flagHumanBase = 1;
                }
            }

            if(flagHumanBase == 0) continue; // ヒトの塩基がギャップの場合は計算(考慮)しない
            if(flagTestis == 0 && flagEye == 0 && flagBrain == 0) continue;



            /*
                Q[n][m][m]     => Q[i][j][k] : 頂点iから親方向への区間(辺)上で塩基置換が(1回以上)起こり、親で塩基k、時刻t_iにおいて塩基jになっている確率。
                                            ※行列計算のため計算は別で行う(scipy.linalg.eig)。

                alpha[n][m] => alpha[i][j] : 頂点i(塩基j)において(子方向)部分木が整合性(与えられる塩基条件)を満たしている確率。
                                            つまり、葉からボトムアップにDPをする。
                                            よって、 
                                                alpha[i][j] += (alpha[V[i].sonLeft][k]*(塩基jが間隔V[i].distLeftで塩基kになる確率)) * (alpha[V[i].sonRight][l]*(塩基jが間隔V[i].distRightで塩基lになる確率))  (for k,lの二重ループ)
                                            となる。                          
                                            
                beta[n][m]  => beta[i][j]  : あるヒトのパス上の頂点i(塩基j)において、ヒトのパスじゃないほうの子部分が整合性を満たし、かつ親方向の残りの部分も整合性を満たす確率。
                                            つまり、alphaを先に計算してから、ヒトのパス上の頂点を根から葉方向に辿ってDPする。(ヒトのパス上以外の点は計算しない。)
                                            よって、
                                                beta[i][j] += (alpha[ヒトじゃないほうの子][k]*(塩基jがiと子の間で塩基kになる確率)) * (beta[V[i].parent][l]*(塩基lが間隔V[i].distParentで塩基jになる確率))  (for k,lの二重ループ)
                                            となる。

                gamma[n][m] => gamma[i][j] : あるヒトのパス上の頂点i(塩基j)の部分木において、ヒトのパス上ではもう塩基置換は起こらず、かつ部分木全体の整合性を満たす確率。
                                            つまり、alpha, betaを先に計算してから、ヒトのパス上の頂点を葉から根方向に辿ってDPする。(ヒトのパス上以外の点は計算しない。)
                                            よって、
                                                gamma[i][j] += (gamma[ヒトパス上の子][j]*(塩基jがiと子の間で置換しない確率)) * (alpha[ヒトじゃないほうの子][k]*(塩基jがiと子の間で塩基kになる確率))  (for k = 0..(m-1))
                                            となる。

                double Z                   : 系統樹および初期条件qが与えられたもとでの葉の塩基アライメントカラムXの出力確率。
                                            つまり、q * alpha[root]。

                nu[n][m]    => nu[i][j]    : 葉頂点iにおいて、jがアライメントカラムの塩基(X[i]、V[i].base)と一致していたら1、そうでなければ0。ただし、アライメントカラムがgap(-)の場合1。
                                            つまり、V[i].base = 0(塩基はA)のときnu[i] = {1,0,0,0}、v[i].base = 1(塩基はG)のときnu[i] = {0,1,0,0}、　v[i].base = -1(gap)のときnu[i] = {1,1,1,1}などとなる。

                prob[n][m]  => prob[i][j]  : 求めたい確率。条件Zにおいて、あるヒトのパス上の頂点i(塩基j)から親方向への区間(辺)上でMRS(the most recent substitution: 最近塩基置換)が起こる条件付き確率。
                                            言い換えると、条件Zにおいてヒトの進化のTMRS(the time to the most recent substitution event : 最近塩基置換時刻)が辺i上にある確率。
                                            つまり、
                                                prob[i][j] += gamma[i][j] * Q[i][k] * beta[V[i].parent][k]  (for k = 0..(m-1))
                                            となる。

                Pt[n][m][m]                : 行列 U * exp(t*lambda) * U_inv (tは枝長)
                                            自由に塩基置換が起こる範囲
            */
            vector<vector<vector<double>>> Q(n, vector<vector<double>>(m, vector<double>(m)));
            vector<vector<double>> alpha(n, vector<double>(m));
            vector<vector<double>> beta(n, vector<double>(m));
            vector<vector<double>> gamma(n, vector<double>(m));
            double Z = 0;
            vector<vector<double>> nu(n, vector<double>(m));
            vector<double> prob(n);
            vector<vector<vector<double>>> Pt(n, vector<vector<double>>(m, vector<double>(m)));

            // nuの前計算
            nu = nuInIt(V);

            // Ptの前計算
            rep(i,n) {
                if(V[i].flagRoot) continue;
                Pt[i] = calcPt(V[i].distParent);
            }

            // 動的計画法

            /*
                alpha メモ
                - 葉以外は必ず子を二つ持っている。
            */

            alpha = alphaCalc(V, nu, Pt);

            /*
                beta メモ
                根からDP
            */

            beta = betaCalc(V, Pt, alpha);

            /*
                gamma メモ
                humanからDP
            */
            gamma = gammaCalc(V, nu, Pt, alpha);
            


            // Zの計算
            Z = calcZ(V, alpha);

            /*
                Q メモ
            */
            rep(i,n) {
                if(!V[i].flagHumanPath || V[i].flagRoot) continue;
                Q[i] = calcQ(V[i]);
            }   


            /*
                prob メモ
            */
            prob = calcProb(V, Q, beta, gamma, Z);

            vector<double> probabilities;
            double probSum = 0;
            double ma = -1;
            int argma = -1;
            rep(i,n) {
                if(!V[i].flagHumanPath) continue;
                probSum += prob[i];
                probabilities.push_back(prob[i]);
                if(ma < prob[i]) {
                    ma = prob[i];
                    argma = i;
                }
            }

            
                

            int gapCut = 1;
            if(gapCut){
                if(group[maxUnGap] >= argma) continue;
            }


            if(flagTestis && loop % 20 == 0){
                string cc;
                rep(i,numOfLeaf) {
                    char c = enhanserSequences[i][loop];
                    c = toupper(c);
                    // if(c == '-' || c == '.') continue;
                    cout << c;
                    cc += c;
                }
                cout << ' ';
                cout << ' ' << argma << ' ' << group[maxUnGap] << ' ' << ma;
                // cout << ' ' << argma << ' ' << ma << endl;
                // int stoop; cin >> stoop;
                char cNow = '*';
                int countChangeCharacter = 0;
                rep(i,cc.size()){
                    if(cc[i] != cNow){
                        cNow = cc[i];
                        countChangeCharacter++;
                    }
                }
                cout << " Testis " << countChangeCharacter << endl;
                countTestis++;
                averageTestis += countChangeCharacter;
            }
            /*
            if(flagEye){
                string cc;
                rep(i,numOfLeaf) {
                    char c = enhanserSequences[i][loop];
                    c = toupper(c);
                    if(c == '-' || c == '.') continue;
                    cout << c;
                    cc += c;
                }
                cout << ' ';
                cout << ' ' << argma << ' ' << group[maxUnGap] << ' ' << ma;
                // cout << ' ' << argma << ' ' << ma << endl;
                // int stoop; cin >> stoop;
                char cNow = '*';
                int countChangeCharacter = 0;
                rep(i,cc.size()){
                    if(cc[i] != cNow){
                        cNow = cc[i];
                        countChangeCharacter++;
                    }
                }
                cout << " Eye " << countChangeCharacter << endl;
                countEye++;
                averageEye += countChangeCharacter;
            }
            if(flagBrain){
                string cc;
                rep(i,numOfLeaf) {
                    char c = enhanserSequences[i][loop];
                    c = toupper(c);
                    if(c == '-' || c == '.') continue;
                    cout << c;
                    cc += c;
                }
                cout << ' ';
                cout << ' ' << argma << ' ' << group[maxUnGap] << ' ' << ma;
                // cout << ' ' << argma << ' ' << ma << endl;
                // int stoop; cin >> stoop;
                char cNow = '*';
                int countChangeCharacter = 0;
                rep(i,cc.size()){
                    if(cc[i] != cNow){
                        cNow = cc[i];
                        countChangeCharacter++;
                    }
                }
                cout << " Brain " << countChangeCharacter << endl;
                countBrain++;
                averageBrain += countChangeCharacter;
            }
            */

            argmaxProb.push_back(argma);
            maxProbScore.push_back(ma);
            allProbabilities.push_back(probabilities);
        }

        // ファイル書き込み
        writeToOfstream(enhancerName, tissueCount, tissues, allProbabilities, ofs);

        // 標準出力
        /*
        cout << enhancerIte << endl;
        cout << enhancerName << endl;
        cout << tissueCount << endl;
        rep(i,tissueCount) cout << tissues[i] << endl;

        int countHumanPath = 19; // 数えない
        int countArgMaxProb[countHumanPath] = {};
        rep(i,argmaxProb.size()) {
            countArgMaxProb[argmaxProb[i]]++;
        }
        rep(i,countHumanPath) cout << setfill('0') << right << setw(3) << countArgMaxProb[i] << ' ';
        cout << endl;
        cout << endl;
        */
        

        enhancerIte++;
    }

    averageTestis /= countTestis;
    averageEye /= countEye;
    averageBrain /= countBrain;

    cout << countEye << ' ' << countTestis << ' ' << countBrain << endl;
    cout << "averageEye = " << averageEye << ", averageTestis = " << averageTestis << ", averageBrain = " << averageBrain << endl;

    end = clock();
    cout << "Total time = " << (double)(end-start)/CLOCKS_PER_SEC << " sec" << endl;
    return 0;
}




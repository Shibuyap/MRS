#include <bits/stdc++.h>
#define rep(i,n) for(int i = 0; i < (n); ++i)
#define srep(i,s,t) for (int i = s; i < t; ++i)
#define drep(i,n) for(int i = (n)-1; i >= 0; --i)
using namespace std;
typedef long long int ll;
typedef pair<int,int> P;
#define yn {puts("Yes");}else{puts("No");}
#define MAX_N 200005

int main() {
    clock_t start,end;
	start = clock();

    // enhancer_tissue.txtのデータは一旦mapに保管する(string型 7*10^3程度)
    string fileNameIfsTissue = "enhancer_tissue.txt";
    const char* cstrIfsTissue = fileNameIfsTissue.c_str();
    ifstream ifsTissue(cstrIfsTissue);
    /*
        1行目はヘッダー
        2~7104行目がエンハンサー情報
        7105行目は改行
        各行は chr~ UBERON~ float型の値 tissue となっている
        各行stringで4つ読み込んで1番目と4番目を使う

        mapはキーがchr~、値は通し番号とし、vector<string> tissues[10000]で器官名を保管
    */
    map<string, int> chrToNum;
    set<string> allTissue; // 器官名列挙用
    vector<string> tissues[10000];
    int chrCount = 0;
    string lineTissue;
    int lineCountTissue = 0;
    while(getline(ifsTissue, lineTissue)){
        lineCountTissue++;
        // if(lineCountTissue == 1) continue;
        rep(_,7103){
            int spaceCount = 0;
            string chr, uberon, val, tissueName;
            ifsTissue >> chr >> uberon >> val >> tissueName;
            // if(_ == 0) cout << chr << ' ' << tissueName << endl;
            /*
            rep(i,lineTissue.length()){
                if(lineTissue[i] == ' '){
                    spaceCount++;
                }else{
                    if(spaceCount == 0) chr += lineTissue[i];
                    if(spaceCount == 1) uberon += lineTissue[i];
                    if(spaceCount == 2) val += lineTissue[i];
                    if(spaceCount == 3) tissueName += lineTissue[i];
                }
            }
            */
            // cout << chr << endl;
            if(chrToNum.find(chr) == chrToNum.end()){
                chrToNum[chr] = chrCount;
                tissues[chrCount].push_back(tissueName);
                chrCount++;
            }else{
                tissues[chrToNum[chr]].push_back(tissueName);
            }
            if(tissueName.size() == 0) {
                cout << _ + 1 << endl;
            }
            allTissue.insert(tissueName);
        }
        break;
    }

    cout << allTissue.size() << endl;
    for(auto t : allTissue) cout << t << endl;
    int stoppp; cin >> stoppp;
    
    
    // mafファイル
    /*
        ## mafファイルについて
        - 全部で3339987行(終わりの改行は消した、本当は3339988行)
        - ##maf version=1 scoring=zero が1行目および途中にあり計100行ある。
        - エンハンサーは各102行で32744個ある
        - 1行目はa score=0.000000、2~101行目は s chr1:858256-858648 0 392 + 392 配列、102行目は改行でワンセット
        - 32744 * 102 + 100 = 3339988となり計算が合う

        ## データ成型
        - mafファイルを1000エンハンサーごとにtxtファイルに分割する。(計33ファイル)
        - 1行目にエンハンサー名(例 chr1:858256-858648)
        - 2行目に対応する器官の数k(kは0以上の整数)
        - 対応する器官が1つ以上あれば3~(2+k)行目まで器官の名称(例 brain) (つまり1器官毎に改行する。)
        - 3+k行目から100行(102+k行目まで)各配列
        - 103+k行目 改行
        - 以上1000セット(33番目のファイルは744セット)
    */
    string fileNameIfsMaf = "enhancers.maf";
    const char* cstrIfsMaf = fileNameIfsMaf.c_str();
    ifstream ifsMaf(cstrIfsMaf);

    int allUseLineOfTissueTxtCount = 0;
    /*
        動かしてみた結果
        allUseLineOfTissueTxtCount = 7102
        Total time = 44.472 sec
        7102が7103と微妙に合わない。配列のないエンハンサー混ざってる？ <- バグな気がしてきた
        後で調べる
    */
    int fileNumber = 0;
    while(true){
        fileNumber++;

        
        
        // enhancers_out_01.txtのようなファイルができる
        string fileNameOut = "enhancers_out_";
        if(fileNumber <= 9){
            fileNameOut += '0';
            fileNameOut += to_string(fileNumber);
        }else{
            fileNameOut += to_string(fileNumber);
        }
        fileNameOut += ".txt";

        const char* cstrOfs = fileNameOut.c_str();
        ofstream ofs(cstrOfs); // ファイル書き込み

        int speciesNumber = 100; // 生物種の数
        int enhancerCount = 0; // 1000まで数える
        string lineMaf; // 1行読み取りに使う

        // デバッグ用
        int cnt = 0;
        int cnts = 0;
        int cnta = 0;
        int cntsharp = 0;

        while (getline(ifsMaf, lineMaf)){
            cnt++;
            if(lineMaf[0] == 's') cnts++;
            if(lineMaf[0] == 'a') cnta++;
            if(lineMaf[0] == '#') cntsharp++;
            
            if(lineMaf[0] == '#'){
                continue;
            } else if(lineMaf[0] == 'a'){ // 次の100行を使う
                rep(i,speciesNumber){
                    cnt++;
                    string t[7]; // 各行空白区切りで7文字列あるため
                    rep(j,7){
                        ifsMaf >> t[j]; 
                    }

                    // ヒトのときだけtissue.txtのデータと照合する
                    // t[1]がchr~
                    if(i == 0){
                        ofs << t[1] << endl;
                        if(chrToNum.find(t[1]) == chrToNum.end()){
                            ofs << 0 << endl;
                        }else{
                            // cout << enhancerCount << ' ' << t[1] << endl;
                            allUseLineOfTissueTxtCount += tissues[chrToNum[t[1]]].size();
                            ofs << tissues[chrToNum[t[1]]].size() << endl;
                            rep(j,tissues[chrToNum[t[1]]].size()){
                                ofs << tissues[chrToNum[t[1]]][j] << endl;
                            }
                        }
                    }

                    // ファイル書き込み
                    ofs << t[6] << endl;
                }
                
                enhancerCount++;
                if(enhancerCount == 1000) break;
            }           
        }

        cout << "fileNumber = " << fileNumber << endl;
        cout << "enhancerCount = " << enhancerCount << endl;
        cout << "cnt = " << cnt << endl;
        cout << "cnts = " << cnts << endl;
        cout << "cnta = " << cnta << endl;
        cout << "cntsharp = " << cntsharp << endl;
        cout << endl;

        ofs.close();

        if(enhancerCount == 0) break;
    }

        
    cout << "OK" << endl;
    cout << "allUseLineOfTissueTxtCount = " << allUseLineOfTissueTxtCount << endl;
    end = clock();
    cout << "Total time = " << (long double)(end-start)/CLOCKS_PER_SEC << " sec" << endl;
    return 0;
}
 
 

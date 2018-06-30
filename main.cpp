#include<string>
#include<ctime>
#include<cmath>

#include<thread>
#include<mutex>
 
#include <map>
#include <fstream>
#include <sstream>
#include <cstdlib>


 const double INF = 1e9;
 const double EPS = 1e-9;
 const double PI = 3.1415926535897932384626;
	const int KD_MAX_THREADS = 5;
	const int PM_MAX_THREADS = 9;

	
	
	
	






extern const int KD_MAX_THREADS;


static int mtxTrees;
std::mutex  *mtx;
class TriangleTree;

typedef unsigned char byte;
typedef unsigned short word;
typedef unsigned int dword;



extern const double INF;
extern const double EPS;
extern const double PI;




	static unsigned mSeed=0;
    inline   double ran() {//RandomLCG,产生0-1之间的随机数 
        mSeed = 214013 * mSeed + 2531011;
		return  mSeed * (1.0 / 4294967296);//快return ((double)mSeed  )/ 4294967296 ;//	return  mSeed * (( double)1 / 4294967296);
	}
	inline   int ran1() {//RandomLCG，产生0-40亿之间的随机数 
        mSeed = 214013 * mSeed + 2531011;
		return mSeed;
	}  
    
    
class ColorRGB {
public:
	float r , g , b;
	ColorRGB( double R = 0 , double G = 0 , double B = 0 ) : r( R ) , g( G ) , b( B ) {}
	~ColorRGB() {}

	friend ColorRGB operator + ( const ColorRGB&A , const ColorRGB&B ){
		return ColorRGB( A.r + B.r , A.g + B.g , A.b + B.b );
	}
	friend ColorRGB operator - ( const ColorRGB&A , const ColorRGB&B ){
		return ColorRGB( A.r - B.r , A.g - B.g , A.b - B.b );
	}
	friend ColorRGB operator * ( const ColorRGB&A , const ColorRGB&B ) {
		return ColorRGB( A.r * B.r , A.g * B.g , A.b * B.b );
	}
	friend ColorRGB operator * ( const ColorRGB&A , const double&k ) {
		return ColorRGB( A.r * k , A.g * k , A.b * k );
	}
	friend ColorRGB operator / ( const ColorRGB&A , const double&k ) {
		return ColorRGB( A.r / k , A.g / k , A.b / k );
	}
	friend ColorRGB& operator += ( ColorRGB&A , const ColorRGB&B ){
		A = A + B;
		return A;
	}
	friend ColorRGB& operator -= ( ColorRGB&A , const ColorRGB&B ){
		A = A - B;
		return A;
	}
	friend ColorRGB& operator *= ( ColorRGB&A , const ColorRGB&B ){
		A = A * B;
		return A;
	}
	friend ColorRGB& operator *= ( ColorRGB&A , const double&k ){
		A = A * k;
		return A;
	}
	friend ColorRGB& operator /= ( ColorRGB&A , const double&k ){
		A = A / k;
		return A;
	}
	
	inline double Average3(){//3个分量r,g,b平均 
		return ( r + g + b ) / 3;
	}
	inline double RGBMax(){//3个分量中最大的一个 
		if (r > g)
			return (r > b) ? r : b;
		return (g > b) ? g : b;
	}
	inline bool IsZero(){//3个分量 r,g,b分别趋近于0 
		return fabs( r ) < EPS && fabs( g ) < EPS && fabs( b ) < EPS;
	}

	inline ColorRGB AmLimit(){ // 3个分量都不大于1；	
		return ColorRGB( std::min( r , 1.0f ) , std::min( g , 1.0f ) , std::min( b , 1.0f ) );
	}
	inline ColorRGB Exp(){//3个分量，分别求指数，e的r次方,e的g次方,e的b次方。	
		return ColorRGB(exp(r), exp(g), exp(b));
	}
	
};

    
class VDouble3 {
public:
	double x , y , z;
	
	VDouble3( double X = 0 , double Y = 0 , double Z = 0 ) : x( X ) , y( Y ) , z( Z ) {}
	~VDouble3() {}

	friend VDouble3 operator + ( const VDouble3&A , const VDouble3&B ){
		return VDouble3( A.x + B.x , A.y + B.y , A.z + B.z );
	}
	friend VDouble3 operator - ( const VDouble3&A , const VDouble3&B ) {
		return VDouble3( A.x - B.x , A.y - B.y , A.z - B.z );
	}
	friend VDouble3 operator * ( const VDouble3&A , const double&k ) {
		return VDouble3( A.x * k , A.y * k , A.z * k );
	}
	friend VDouble3 operator * ( const VDouble3&A , const VDouble3&B ) {
		return VDouble3( A.x * B.x , A.y * B.y , A.z * B.z );
	}
	friend VDouble3 operator / ( const VDouble3&A , const double&k ) {
		return VDouble3( A.x / k , A.y / k , A.z / k );
	}

	friend VDouble3& operator += ( VDouble3&A , const VDouble3&B ){
		A = A + B;
		return A;
	}
	friend VDouble3& operator -= ( VDouble3&A , const VDouble3&B ){
		A = A - B;
		return A;
	}
	friend VDouble3& operator *= ( VDouble3&A , const double&k ) {
		A = A * k;
		return A;
	}
	friend VDouble3& operator *= ( VDouble3& A , const VDouble3& B ) {
		A = A * B;
		return A;
	}
	friend VDouble3& operator /= ( VDouble3& A , const double& k ) {
		A = A / k;
		return A;
	}
	
	
	
	
	
	friend VDouble3 operator - ( const VDouble3& A ) {
		return VDouble3( -A.x , -A.y , -A.z );
	}
	 
	inline double& GetCoord( int axis ){
		if ( axis == 0 ) return x;
		if ( axis == 1 ) return y;
		if ( axis == 2 ) return z;
	}//得到某一个坐标（3坐标中的1个）的值， 
	inline VDouble3 GetUnitVect(){
		double module = Module();
		if (module < EPS) return VDouble3(0, 0, 1);
		if ( fabs(module -1)< EPS) return *this;
		return *this / module;
	}
	inline VDouble3 GetVerticalVect() {//返回 与Z轴叉积 的单位向量 
	//	return  Cross(VDouble3( 0 , 0 , 1 )).GetUnitVect();
		VDouble3 ret = Cross(VDouble3( 0 , 0 , 1 ));
		if ( ret.IsZero() ) return  VDouble3( 1 , 0 , 0 );	
		return ret.GetUnitVect();
	}
	
	 
	inline double Dot( const VDouble3&term ){
		return x * term.x + y * term.y + z * term.z;
	}
	inline VDouble3 Cross(const VDouble3&term) {
		return VDouble3(y * term.z - z * term.y , z * term.x - x * term.z , x * term.y - y * term.x );
	}
	
	
	inline double Module2() {
		return x * x + y * y + z * z;
	}
	inline double Module(){
		return sqrt( x * x + y * y + z * z );
	}
	inline double Distance2( VDouble3&term){
		return ( term - *this ).Module2();
	}
	inline double Distance( VDouble3&term ) {
		return ( term - *this ).Module();
	}
	
	
	inline bool IsZero(){
		return fabs( x ) < EPS && fabs( y ) < EPS && fabs( z ) < EPS;
	}
	inline void AssRandomVect() {
		double locadist2;
		do 	{
			x = 2 * ran() - 1;
			y = 2 * ran() - 1;
			z = 2 * ran() - 1;
			locadist2=x * x + y * y + z * z;
		} while (locadist2> 1 || locadist2< EPS );
		locadist2=sqrt(locadist2);
		x/= locadist2;
		y/= locadist2;
		z/= locadist2;
	}
	inline void Input( std::stringstream& fin ) {
		fin >> x >> y >> z;
	}
	
	
	inline VDouble3 Diffuse() {//theta   phi
		return Rotate( GetVerticalVect(), acos( sqrt( ran() ) ) ).Rotate( *this , ran() * 2 * PI );
	}
	inline VDouble3 Reflect( VDouble3 N ) {
		return *this - N * ( 2 * Dot( N ) );  //计算该光线遇到面N的反射光线
	}
	inline VDouble3 Refract( VDouble3 N , double zhe , bool* refracted ) {
		VDouble3 In = GetUnitVect();
		double cosA1 = -N.Dot( In ) ;
		double cosA2 = 1 -  zhe * zhe * ( 1 - cosA1 * cosA1 );
		if ( cosA2 <= EPS ) return In.Reflect( N );
		if (refracted != NULL) *refracted ^= true;
		return N * ( zhe * cosA1 - sqrt( cosA2 ) ) +In * zhe ;	
	}
	inline VDouble3 Rotate( VDouble3 axis , double theta ) {
		VDouble3 ret;
		double cost = cos( theta );
		double sint = sin( theta );
		ret.x += x * ( axis.x * axis.x + ( 1 - axis.x * axis.x ) * cost );
		ret.x += y * ( axis.x * axis.y * ( 1 - cost ) - axis.z * sint );
		ret.x += z * ( axis.x * axis.z * ( 1 - cost ) + axis.y * sint );
		ret.y += x * ( axis.y * axis.x * ( 1 - cost ) + axis.z * sint );
		ret.y += y * ( axis.y * axis.y + ( 1 - axis.y * axis.y ) * cost );
		ret.y += z * ( axis.y * axis.z * ( 1 - cost ) - axis.x * sint );
		ret.z += x * ( axis.z * axis.x * ( 1 - cost ) - axis.y * sint );
		ret.z += y * ( axis.z * axis.y * ( 1 - cost ) + axis.x * sint );
		ret.z += z * ( axis.z * axis.z + ( 1 - axis.z * axis.z ) * cost );
		return ret;
	}
};




struct BITMAPFILEHEADER {
	dword bfSize;
	word bfReserved1;
	word bfReserved2;
	dword bfOffBits;
};

struct BITMAPINFOHEADER {
	dword biSize;
	long biWidth;
	long biHeight;
	word biPlanes;
	word biBitCount;
	dword biCompression;
	dword biSizeImage;
	long biXPelsPerMeter;
	long biYPelsPerMeter;
	dword biClrUsed;
	dword biClrImportant;
};

struct RGBQUAD {
	byte rgbBlue;
	byte rgbGreen;
	byte rgbRed;
	byte rgbReserved;
};

struct IMAGEDATA {
	byte red;
	byte green;
	byte blue;
	ColorRGB GetColor255() {
		return ColorRGB( red/ 255.0 , green/ 255.0 , blue/ 255.0 ) ;
	}
};

class Bitmap {
	BITMAPFILEHEADER bit_Head;
	BITMAPINFOHEADER bit_Info;
	IMAGEDATA** ima;
public:
	Bitmap( int H = 0 , int W = 0 )  {
		Initialize( H , W );
	}
	~Bitmap(){
		if(ima==NULL)return;
		for ( int i = 0 ; i < bit_Info.biHeight ; i++ ) delete[] ima[i];
		delete[] ima;
		ima=NULL;
	}


	void Initialize( int H , int W ) {
		bit_Head.bfReserved1 = 0;
		bit_Head.bfReserved2 = 0;
		bit_Head.bfOffBits = 54;
	
		bit_Info.biSize = 40;
		bit_Info.biPlanes = 1;
		bit_Info.biHeight = H;
		bit_Info.biWidth = W;
		bit_Info.biBitCount = 24;
		bit_Info.biCompression = 0;
		bit_Info.biSizeImage = H * W * 3;
		bit_Info.biXPelsPerMeter = 0;
		bit_Info.biYPelsPerMeter = 0;
		bit_Info.biClrUsed = 0;
		bit_Info.biClrImportant = 0;
	
		bit_Head.bfSize = bit_Info.biSizeImage + bit_Info.biBitCount;
		if(W<=0||H<=0){
			ima=NULL;
			return;
		}
		ima = new IMAGEDATA*[H];
		for ( int i = 0 ; i < H ; i++ )
			ima[i] = new IMAGEDATA[W];
		printf("head54=%d,w64=%d,h48=%d",sizeof(bit_Head),W,H);
	}
	void Input( std::string file ) {
		if(ima!=NULL)
		{	for ( int i = 0 ; i < bit_Info.biHeight ; i++ ) delete[] ima[i];
			delete[] ima;
			ima=NULL;
		}
		FILE *fpi = fopen( file.c_str() , "rb" );
		if(fpi==NULL)return;
		word bfType;
		fread( &bfType , 1 , 2 , fpi );
		fread( &bit_Head , 1 , sizeof( BITMAPFILEHEADER ) , fpi );
		fread( &bit_Info , 1 , sizeof( BITMAPINFOHEADER ) , fpi );
		
	
		Initialize( bit_Info.biHeight , bit_Info.biWidth );	
		int number = bit_Info.biWidth % 4;
		for(int i = 0 ; i < bit_Info.biHeight ; i++ ) {
			for(int j = 0 ; j < bit_Info.biWidth ; j++ ) {
				fread( &ima[i][j].blue , 1 , 1 , fpi );
				fread( &ima[i][j].green , 1 , 1 , fpi );
				fread( &ima[i][j].red , 1 , 1 , fpi );
			}
			byte buffer = 0;
			for (int j = 0; j < number; j++) fread( &buffer , 1 , 1 , fpi );
		}
		fclose( fpi );
	}

	void Output( std::string file ) {
		printf("bmp1  out head54=%d, w=%d,h=%d\n", sizeof( BITMAPFILEHEADER )+sizeof( BITMAPINFOHEADER )+2
			,bit_Info.biWidth,  bit_Info.biHeight );
		//	printf("head54=%d,w64=%d,h48=%d",sizeof(BITMAPFILEHEADER),W,H);
		FILE *fpw = fopen( file.c_str() , "wb" );
		if(fpw==NULL)return;
		word bfType = 0x4d42;
		fwrite( &bfType , 1 , 2 , fpw );
		fwrite( &bit_Head , 1 , sizeof( BITMAPFILEHEADER ) , fpw );
		fwrite( &bit_Info , 1 , sizeof( BITMAPINFOHEADER ) , fpw );	
		int number = bit_Info.biWidth % 4;
	
	
		printf("bmp2  out head54=%d, w=%d,h=%d   number=%d\n", sizeof( BITMAPFILEHEADER )+sizeof( BITMAPINFOHEADER )+2
			,bit_Info.biWidth,  bit_Info.biHeight   ,number );
		
	
		for ( int i = 0 ; i < bit_Info.biHeight ; i++ ) {
			for ( int j = 0 ; j < bit_Info.biWidth ; j++ ) {
				fwrite( &ima[i][j].blue , 1 , 1 , fpw );
				fwrite( &ima[i][j].green , 1 , 1 , fpw );
				fwrite( &ima[i][j].red , 1 , 1 , fpw );
			}
			byte buffer = 0;
			for (int j = 0; j < number; j++) 
			{
				printf("bmp3333 number=%d\n",  number );
	
				fwrite( &buffer , 1 , 1 , fpw );
			}
		}
		fclose( fpw );
	}
	
	inline int GetH() { return bit_Info.biHeight; }
	inline int GetW() { return bit_Info.biWidth; }

	inline void SetColor( int i , int j , ColorRGB col ){
		ima[i][j].red = ( int ) ( col.r * 255 );
		ima[i][j].green = ( int ) ( col.g * 255 );
		ima[i][j].blue = ( int ) ( col.b * 255 );
	}
	ColorRGB GetSmoothColor( double u , double v ) {
		double U = ( u + EPS - floor( u + EPS ) ) * bit_Info.biHeight;
		double V = ( v + EPS - floor( v + EPS ) ) * bit_Info.biWidth;
		int U1 = ( int ) floor( U + EPS ) , U2 = U1 + 1;
		int V1 = ( int ) floor( V + EPS ) , V2 = V1 + 1;
		double rat_U = U2 - U;
		double rat_V = V2 - V;
		if ( U1 < 0 ) U1 = bit_Info.biHeight - 1; 
		if ( U2 == bit_Info.biHeight ) U2 = 0;
		if ( V1 < 0 ) V1 = bit_Info.biWidth - 1; 
		if ( V2 == bit_Info.biWidth ) V2 = 0;
		ColorRGB ret;
		ret = ret + ima[U1][V1].GetColor255() * rat_U * rat_V;
		ret = ret + ima[U1][V2].GetColor255() * rat_U * ( 1 - rat_V );
		ret = ret + ima[U2][V1].GetColor255() * ( 1 - rat_U ) * rat_V;
		ret = ret + ima[U2][V2].GetColor255() * ( 1 - rat_U ) * ( 1 - rat_V );
		return ret;
	}
	
};




class MaterialQua {
public:
	ColorRGB color , absor;
	double refl , refr;
	double diff , spec;
	double rindex;
	double drefl;
	Bitmap* texture;

	MaterialQua(ColorRGB color1 , ColorRGB absor1,
		double refl1 , double refr1,double diff1 , double spec1,
		double rindex1,double drefl1,std::string texture1 	) {
		color =color1; 
		absor = absor1;//ColorRGB();
		refl = refl1;
		refr = refr1;
		diff = diff1;
		spec = spec1;
		rindex = rindex1;
		drefl = drefl1;
		if(texture1==""){
			texture=NULL;
			return ;
		}
		texture = new Bitmap;
		texture->Input( texture1 );
	 
	
		printf("---w=%d,h=%d\n", texture->GetW(), texture->GetH());
	}
	~MaterialQua() {
		if(texture==NULL)return ;
		delete texture;
		texture=NULL;
	}
	double BRDF(VDouble3 ray_R, VDouble3 N, VDouble3 ray_I) {
		double ret = 0;
		ray_R = ray_R.GetUnitVect();
		ray_I = ray_I.GetUnitVect();
		
		if (diff > EPS && ray_R.Dot(N) > EPS)
			ret += diff * ray_R.Dot(N);
		if (spec > EPS && ray_R.Dot(-ray_I.Reflect(N)) > EPS)
			ret += spec * pow(ray_R.Dot(-ray_I.Reflect(N)), 50);
		return ret;
	}
};


//本文件含有5个类，是最基础的类：光线类Ray，碰撞记录类ColliderRecord， 基础物体类BaseObject， 
//非光源物体类NoLightObj ，光源类LightObject 

class BaseObject;
 
class Ray {	
public:
	VDouble3 dir; //位置和入射方向
	VDouble3 pos; 
	ColorRGB power;//光子能量
	int plane;//光子在KD树中划分的平面
	Ray(VDouble3 pos1,VDouble3 dir1 ):  pos(pos1),dir(dir1) {}
	Ray(VDouble3 pos1  ):  pos(pos1)  {}
	Ray( )  {}
	~Ray() {}
};
class ColliderRecord {
public:
	BaseObject* priBs;//碰撞到的物体
	double dista;//光线到碰撞点的距离---碰撞前光线走过距离
	bool impact, frontage;//impact是否碰撞    frontage是否在物体对象正面碰撞
	VDouble3 No, Co, In;//碰撞点法向量No   碰撞点坐标Co    入射光线的方向In 

	ColliderRecord() {
		priBs = NULL;
		dista=0;
		impact = false;
	}
	~ColliderRecord() {}
};
class BaseObject {
public:
	int hash_d;//一个hash随机数
	BaseObject* next;//下一个指针，这里不用删除内存（用链表储存，在场景中被删除）
	BaseObject(){
		hash_d = ran1();
		next = NULL;
	}
	~BaseObject() {}	
	int getHash_d() { return hash_d; }
	BaseObject* GetNext() { return next; }
	void SetNext( BaseObject* light ) { next = light; }

	virtual ColorRGB GetColor255() { return ColorRGB(); }
	virtual void SetMaterial(MaterialQua* _material) {  }
	virtual MaterialQua* GetMaterial() { return NULL; }

	virtual ColliderRecord CollideOne(  Ray ray ) {return ColliderRecord();}//计算光线到物体对象的碰撞情况，返回碰撞记录 
	//virtual VDouble3 GetAAA() = 0;
};

class NoLightObj: public BaseObject {
public:
	MaterialQua* material;
	NoLightObj(){}
	
	NoLightObj(ColorRGB color1 , ColorRGB absor1,
		double refl1 , double refr1,double diff1 , double spec1,
		double rindex1,double drefl1,std::string texture1 	):BaseObject(){
			
		material = new MaterialQua( color1 ,  absor1,
		 refl1 ,  refr1, diff1 ,  spec1,
		 rindex1, drefl1, texture1 	);
		
	}
	~NoLightObj(){
		delete material;
	}

	void SetMaterial(MaterialQua* _material) { material = _material; }
	MaterialQua* GetMaterial() { return material; }
	
	virtual ColliderRecord CollideOne(  Ray ray ){return ColliderRecord();}
	
	virtual void PreTreatment() {}//在读入完后，进行预处理	
	virtual ColorRGB GetTexture(VDouble3 C){
		return ColorRGB(0, 0, 0);
	}//碰撞点的纹理颜色
};
class LightObject: public BaseObject {
public:
	ColorRGB color;//光源的颜色
	LightObject():BaseObject(){ }
	~LightObject() {}	
	
	ColorRGB GetColor255() { return color; }
	
	virtual ColliderRecord CollideOne(  Ray ray ) = 0;


	virtual ColorRGB getIrra( ColliderRecord* colli_R , NoLightObj* NObj_head , int shade_quality , int* hash ) = 0;
	//碰撞点的光照颜色
	virtual Ray getPhoton() = 0;//getPhoton   得到一个出来的光子
};

class Near_pho {
public:
	VDouble3 pos;//求关于点pos的最近k光子
	int max_photons , found;//最大保留光子数,已经找到的光子数

	bool got_heap;//是否已建堆
	double* dist2;//最近光子与点pos平方距离
	Ray** photons;

	Near_pho() {
		max_photons = found = 0;
		got_heap = false;
		dist2 = NULL;
		photons = NULL;
	}
	~Near_pho() {
		delete[] dist2;
		delete[] photons;
	}

};





class Polyhedron : public NoLightObj {
	VDouble3 O, size, angles;
	VDouble3* vertexN;
	std::pair<double, double>* pixel;
	std::string meshFile;
	TriangleTree* tree;

public:

	Polyhedron(VDouble3 N1 , VDouble3 Dx1 ,VDouble3 Dy1,double R1,ColorRGB color1 , ColorRGB absor1,
		double refl1 , double refr1,double diff1 , double spec1,
		double rindex1,double drefl1,std::string texture1 );
	~Polyhedron();
	
	VDouble3 GetO() { return O; }
	VDouble3 GetSize() { return size; }
	VDouble3 GetAngles() { return angles; }
	void SetVertexN(VDouble3* _vertexN) { vertexN = _vertexN; }
	VDouble3& GetVertexN(int i) { return vertexN[i]; }
	TriangleTree* GetTree() { return tree; }
	void SetPixel(std::pair<double, double>* _pixel) { pixel = _pixel; }
	std::pair<double, double>& GetPixel(int i) { return pixel[i]; }
	
	
	void PreTreatment();
	ColliderRecord CollideOne(Ray ray);
};


class Triangle : public NoLightObj {
public:
	Polyhedron* parent;
	VDouble3 N, pos[3],B,C;
	int vertex[3], textureVertex[3], normalVectorID[3];
	double dot11, dot01, dot00,inverDeno;
	Triangle(){
		parent = NULL;
		vertex[0] = vertex[1] = vertex[2] = 0;
		textureVertex[0] = textureVertex[1] = textureVertex[2] = 0;
		normalVectorID[0] = normalVectorID[1] = normalVectorID[2] = 0;
	}
	~Triangle() {}

	void SetParent(Polyhedron* _parent) { parent = _parent; }
	VDouble3& GetPos(int i) { return pos[i]; }
	int& GetVertex(int i) { return vertex[i]; }
	int& GetTextureVertex(int i) { return textureVertex[i]; }
	int& GetNormalVectorID(int i) { return normalVectorID[i]; }
	VDouble3& GetN() { return N; }
	
	void Input(std::string var, std::stringstream& fin) {
		if (var == "P0=") pos[0].Input(fin);
		if (var == "P1=") pos[1].Input(fin);
		if (var == "P2=") pos[2].Input(fin);
	}
	void PreTreatment() {
		B = pos[2] - pos[0], C = pos[1] - pos[0];
		dot11 = C.Dot( C);  	dot01 = C.Dot( B);
	  	dot00 = B.Dot( B);
	 	inverDeno = dot01 * dot01 -  dot11 * dot00; 
	 	if(fabs( inverDeno )<EPS ) return  ;   	
		N = C.Cross(B);
		if (N.IsZero()) {
			N = VDouble3(0, 0, 1);
			return;
		}
		N = N.GetUnitVect();
	}
	ColliderRecord CollideOne(Ray ray) {
		ColliderRecord colli_R;	
	   	if(fabs( inverDeno )<EPS ) return colli_R ;    	
		ray.dir = ray.dir.GetUnitVect();
		double nomDi = N.Dot( ray.dir );	
		if ( fabs( nomDi ) <= EPS ) return colli_R; 
		double tt = (  pos[0]-ray.pos).Dot( N ) / nomDi;
		if ( tt <= EPS ) return colli_R; 	
		colli_R.Co = ray.pos + ray.dir * tt;   
		VDouble3 v2 = colli_R.Co -  pos[0];    	
	  	double dot12 = v2.Dot( C);
	  	double dot02 = v2.Dot( B);
	  	double x = (dot01 * dot02 -  dot00 * dot12) / inverDeno; 
	  	if (x < EPS )    return  colli_R;
	  	double y = (dot01 * dot12 -  dot11 * dot02) / inverDeno; 
	  	if (y < EPS || (x + y) > 1.0f-EPS)    return  colli_R; 
	 
	 
		if (parent != NULL && !parent->GetVertexN(normalVectorID[0]).IsZero())
			colli_R.No = parent->GetVertexN(normalVectorID[0]) * (1 - x - y) + parent->GetVertexN(normalVectorID[1]) * x + parent->GetVertexN(normalVectorID[2]) * y;
		else
			colli_R.No = N;
		double d = colli_R.No.Dot(ray.dir);
		colli_R.impact = true;
		colli_R.In = ray.dir;
		colli_R.priBs=this;
		colli_R.dista = tt;
		colli_R.frontage = (d < 0);
		if (!colli_R.frontage) colli_R.No = -colli_R.No;
		return colli_R;
	}
	ColorRGB GetTexture(VDouble3 C) {
		if (material->texture != NULL)return ColorRGB(1,1,1);
		double totalArea = (pos[1] - pos[0]).Cross(pos[2] - pos[0]).Module();
		double area1 = (C - pos[0]).Cross(pos[2] - pos[0]).Module();
		double area2 = (C - pos[0]).Cross(pos[1] - pos[0]).Module();
		double x = 0, y = 0;
		if (totalArea != 0) {
			x = area1 / totalArea;
			y = area2 / totalArea;
		}
		std::pair<double, double> p0 = parent->GetPixel(textureVertex[0]);
		std::pair<double, double> p1 = parent->GetPixel(textureVertex[1]);
		std::pair<double, double> p2 = parent->GetPixel(textureVertex[2]);
		p0.first = p0.first + (p1.first - p0.first) * x + (p2.first - p0.first) * y;
		p0.second = p0.second + (p1.second - p0.second) * x + (p2.second - p0.second) * y;
		return material->texture->GetSmoothColor(p0.first, p0.second);
	}
	double GetMinCoord(int coord) {
		double x0 = pos[0].GetCoord(coord);
		double x1 = pos[1].GetCoord(coord);
		double x2 = pos[2].GetCoord(coord);
		if (x0 < x1)
			return (x0 < x2) ? x0 : x2;
		return (x1 < x2) ? x1 : x2;
	}
	double GetMaxCoord(int coord) {
		double x0 = pos[0].GetCoord(coord);
		double x1 = pos[1].GetCoord(coord);
		double x2 = pos[2].GetCoord(coord);
		if (x0 > x1)
			return (x0 > x2) ? x0 : x2;
		return (x1 > x2) ? x1 : x2;
	}
	
	
};


class TriangleBox {
public:
	VDouble3 minPos, maxPos;
	
	TriangleBox(){
		minPos = VDouble3(INF, INF, INF);
		maxPos = VDouble3(-INF, -INF, -INF);
	}
	~TriangleBox() {}
	
	void Update(Triangle* tri) {
		for (int coord = 0; coord < 3; coord++) {
			if (tri->GetMinCoord(coord) < minPos.GetCoord(coord)) minPos.GetCoord(coord) = tri->GetMinCoord(coord);
			if (tri->GetMaxCoord(coord) > maxPos.GetCoord(coord)) maxPos.GetCoord(coord) = tri->GetMaxCoord(coord);
		}
	}
	bool Cantain(VDouble3 O) {
		for (int coord = 0; coord < 3; coord++)
			if (O.GetCoord(coord) <= minPos.GetCoord(coord) - EPS || O.GetCoord(coord) >= maxPos.GetCoord(coord) + EPS) return false;
		return true;
	}
	double CalnArea() {
		double a = maxPos.x - minPos.x;
		double b = maxPos.y - minPos.y;
		double c = maxPos.z - minPos.z;
		return 2 * (a * b + b * c + c * a);
	}
	double CollideOne(Ray ray) {
		double minDist = -EPS,modu=ray.dir.Module();
		for (int coord = 0; coord < 3; coord++) {
			double times = -1;
			if (ray.dir.GetCoord(coord) >= EPS)
				times = (minPos.GetCoord(coord) -ray.pos.GetCoord(coord)) / ray.dir.GetCoord(coord);
			if (ray.dir.GetCoord(coord) <= -EPS)
				times = (maxPos.GetCoord(coord) -ray.pos.GetCoord(coord)) / ray.dir.GetCoord(coord);
			if (times >= EPS) {
				VDouble3 C =ray.pos + ray.dir * times;
				if (Cantain(C)) {
					double dist =modu * times;//ray.pos.Distance(C);// ( term - *this ).Module(); 
					if (minDist <= -EPS || dist < minDist)
						minDist = dist;
				}
			}
		}
		return minDist;
	}

};


class TriangleNode {
public:
	Triangle** tris;
	int size, plane;//size是三角面的个数 
	double split;
	TriangleBox box;
	TriangleNode* leftNode;
	TriangleNode* rightNode;

	TriangleNode(){
		size = 0;
		plane = -1;
		split = 0;
		leftNode = rightNode = NULL;
	}
	~TriangleNode() {
		for (int i = 0; i < size; i++)
			delete tris[i];
		delete tris;
		delete leftNode;
		delete rightNode;
	}
};


class TriangleTree {
public:
	TriangleNode* root;
	TriangleTree() {
		root = new TriangleNode;
	}
	~TriangleTree() {
		DeleteTree(root);
	}
	TriangleNode* GetRoot() { return root; }

	void DeleteTree(TriangleNode* node) {
		if (node->leftNode != NULL)
			DeleteTree(node->leftNode);
		if (node->rightNode != NULL)
			DeleteTree(node->rightNode);
		delete node;
	}
	void SortTriangle(Triangle** tris, int l, int r, int coord, bool minCoord) {
		double (Triangle::*GetCoord)(int) = minCoord ? &Triangle::GetMinCoord : &Triangle::GetMaxCoord;
		if (l >= r) return;
		int i = l, j = r;
		Triangle* key = tris[(l + r) >> 1];
		while (i <= j) {
			while (j >= l && (key->*GetCoord)(coord) < (tris[j]->*GetCoord)(coord)) j--;
			while (i <= r && (tris[i]->*GetCoord)(coord) < (key->*GetCoord)(coord)) i++;
			if (i <= j) {
				std::swap(tris[i], tris[j]);
				i++;
				j--;
			}
		}
		SortTriangle(tris, i, r, coord, minCoord);
		SortTriangle(tris, l, j, coord, minCoord);
	}
	void DivideNode(TriangleNode* node) {
		if (node->size * KD_MAX_THREADS >= root->size) {
			printf("KDtree(size = %d)\n", node->size);
			mtx->lock();
			mtxTrees++;
			mtx->unlock();
		}
		//iff area0 * size0 + area1 * size1 + totalArea <= totalArea * totalSize then divide
		Triangle** minNode = new Triangle*[node->size];
		Triangle** maxNode = new Triangle*[node->size];
		for (int i = 0; i < node->size; i++) {
			minNode[i] = node->tris[i];
			maxNode[i] = node->tris[i];
		}
		
		double thisCost = node->box.CalnArea() * (node->size - 1);
		double minCost = thisCost;
		int bestCoord = -1, leftSize = 0, rightSize = 0;
		double bestSplit = 0;
		for (int coord = 0; coord < 3; coord++) {
			SortTriangle(minNode, 0, node->size - 1, coord, true);
			SortTriangle(maxNode, 0, node->size - 1, coord, false);
			TriangleBox leftBox = node->box;
			TriangleBox rightBox = node->box;
	
			int j = 0;
			for (int i = 0; i < node->size; i++) {
				double split = minNode[i]->GetMinCoord(coord);
				leftBox.maxPos.GetCoord(coord) = split;
				rightBox.minPos.GetCoord(coord) = split;
				for ( ; j < node->size && maxNode[j]->GetMaxCoord(coord) <= split + EPS; j++);
				double cost = leftBox.CalnArea() * i + rightBox.CalnArea() * (node->size - j);
				if (cost < minCost) {
					minCost = cost;
					bestCoord = coord;
					bestSplit = split;
					leftSize = i;
					rightSize = node->size - j;
				}
			}
	
			j = 0;
			for (int i = 0; i < node->size; i++) {
				double split = maxNode[i]->GetMaxCoord(coord);
				leftBox.maxPos.GetCoord(coord) = split;
				rightBox.minPos.GetCoord(coord) = split;
				for ( ; j < node->size && minNode[j]->GetMinCoord(coord) <= split - EPS; j++);
				double cost = leftBox.CalnArea() * j + rightBox.CalnArea() * (node->size - i);
				if (cost < minCost) {
					minCost = cost;
					bestCoord = coord;
					bestSplit = split;
					leftSize = j;
					rightSize = node->size - i;
				}
			}
		}
	
		delete minNode;
		delete maxNode;
	
		if (bestCoord != -1) {
			leftSize = rightSize = 0;
			for (int i = 0; i < node->size; i++) {
				if (node->tris[i]->GetMinCoord(bestCoord) <= bestSplit - EPS || node->tris[i]->GetMaxCoord(bestCoord) <= bestSplit + EPS)
					leftSize++;
				if (node->tris[i]->GetMaxCoord(bestCoord) >= bestSplit + EPS || node->tris[i]->GetMinCoord(bestCoord) >= bestSplit - EPS)
					rightSize++;
			}
			TriangleBox leftBox = node->box;
			TriangleBox rightBox = node->box;
			leftBox.maxPos.GetCoord(bestCoord) = bestSplit;
			rightBox.minPos.GetCoord(bestCoord) = bestSplit;
			double cost = leftBox.CalnArea() * leftSize + rightBox.CalnArea() * rightSize;
	
			if (cost < thisCost) {
				node->plane = bestCoord;
				node->split = bestSplit;
	
				node->leftNode = new TriangleNode;
				node->leftNode->box = node->box;
				node->leftNode->box.maxPos.GetCoord(node->plane) = node->split;
				
				node->rightNode = new TriangleNode;
				node->rightNode->box = node->box;
				node->rightNode->box.minPos.GetCoord(node->plane) = node->split;
				
				node->leftNode->tris = new Triangle*[leftSize];
				node->rightNode->tris = new Triangle*[rightSize];
				int leftCnt = 0, rightCnt = 0;
				for (int i = 0; i < node->size; i++) {
					if (node->tris[i]->GetMinCoord(node->plane) <= node->split - EPS || node->tris[i]->GetMaxCoord(node->plane) <= node->split + EPS)
						node->leftNode->tris[leftCnt++] = node->tris[i];
					if (node->tris[i]->GetMaxCoord(node->plane) >= node->split + EPS || node->tris[i]->GetMinCoord(node->plane) >= node->split - EPS)
						node->rightNode->tris[rightCnt++] = node->tris[i];
				}
				node->leftNode->size = leftSize;
				node->rightNode->size = rightSize;
	
				if (node->size * KD_MAX_THREADS >= root->size * 2) {
					//DivideNode( node->leftNode);
					std::thread subThread(&TriangleTree::DivideNode, this, node->leftNode);
					subThread.detach();
				} else
					DivideNode(node->leftNode);
				DivideNode(node->rightNode);
			}
		}
		if (node->size * KD_MAX_THREADS >= root->size) {
			mtx->lock();
			mtxTrees--;
			mtx->unlock();
		}
	}
	ColliderRecord TravelTree(TriangleNode* node, Ray ray) {
		if (!node->box.Cantain(ray.pos) && node->box.CollideOne(ray) <= -EPS)
			return ColliderRecord();
	
		if (node->leftNode == NULL && node->rightNode == NULL) {
			ColliderRecord ret;
			for (int i = 0; i < node->size; i++) {
				ColliderRecord colli_R = node->tris[i]->CollideOne(ray);
				if (colli_R.impact && node->box.Cantain(colli_R.Co) && (!ret.impact || colli_R.dista < ret.dista))
					ret = colli_R;
			}
			return ret;
		}
		
		if (node->leftNode->box.Cantain(ray.pos)) {
			ColliderRecord colli_R = TravelTree(node->leftNode, ray);
			if (colli_R.impact) return colli_R;
			return TravelTree(node->rightNode, ray);
		}
		if (node->rightNode->box.Cantain(ray.pos)) {
			ColliderRecord colli_R = TravelTree(node->rightNode, ray);
			if (colli_R.impact) return colli_R;
			return TravelTree(node->leftNode, ray);
		}
	
		double leftDist = node->leftNode->box.CollideOne(ray);
		double rightDist = node->rightNode->box.CollideOne(ray);
		if (rightDist <= -EPS)
			return TravelTree(node->leftNode, ray);
		if (leftDist <= -EPS)
			return TravelTree(node->rightNode, ray);
		
		if (leftDist < rightDist) {
			ColliderRecord colli_R = TravelTree(node->leftNode, ray);
			if (colli_R.impact) return colli_R;
			return TravelTree(node->rightNode, ray);
		}
		ColliderRecord colli_R = TravelTree(node->rightNode, ray);
		if (colli_R.impact) return colli_R;
		return TravelTree(node->leftNode, ray);
	}
	void BuildTree() {
		mtx = new std::mutex;
		DivideNode(root);
		while (true) {
			mtx->lock();
			if (mtxTrees == 0) break;
			mtx->unlock();
		}
		delete mtx;
	}
	ColliderRecord CollideOne(Ray ray) {
		return TravelTree(root, ray);
	}
	
	
};


class ObjReader {
	Polyhedron* polyhedron;
	int vSize, vtSize, vnSize, fSize, matSize;
	VDouble3* v;
	std::pair<double, double>* vt;
	VDouble3* vn;
	Triangle** tris;
	MaterialQua** mat;
	std::map<std::string, int> matMap;
public:
	ObjReader(){
		polyhedron = NULL;
		vSize = 0;
		vtSize = 0;
		vnSize = 0;
		fSize = 0;
		matSize = 0;
	}
	~ObjReader() {}
	void SetPolyhedron(Polyhedron* _polyhedron) { polyhedron = _polyhedron; }

	
	void ReadMtlSize(std::string file) {
		std::ifstream fin(file.c_str());
		std::string order;
	
		while (getline(fin, order, '\n')) {
			std::stringstream fin2(order);
			std::string var;
			if (!(fin2 >> var)) continue;	
		}
		fin.close();
	
		mat = new MaterialQua*[matSize + 1];
	}	
	void ReadObjSize(std::string file) {
		std::ifstream fin(file.c_str());
		std::string order;
		
		while (getline(fin, order, '\n')) {
			std::stringstream fin2(order);
			std::string var;
			if (!(fin2 >> var)) continue;	
			if (var == "v")
				vSize++;	
			if (var == "f") {
				int vertexCnt = 0;
				std::string var;
				while (fin2 >> var)
					vertexCnt++;
				fSize += std::max(0, vertexCnt - 2);
			}
		}
		fin.close();
	
		v = new VDouble3[vSize + 1];
		vt = new std::pair<double, double>[vtSize + 1];
		if (vnSize == 0)
			vn = new VDouble3[vSize + 1];
		else
			vn = new VDouble3[vnSize + 1];
		tris = new Triangle*[fSize];
	}
	void CalnVn() {
		if (vnSize > 0) {
			for (int i = 1; i <= vnSize; i++) {
				vn[i] = vn[i].Rotate(VDouble3(1, 0, 0), polyhedron->GetAngles().GetCoord(0));
				vn[i] = vn[i].Rotate(VDouble3(0, 1, 0), polyhedron->GetAngles().GetCoord(1));
				vn[i] = vn[i].Rotate(VDouble3(0, 0, 1), polyhedron->GetAngles().GetCoord(2));
			}
		} 
	}	
	void ReadObj(std::string file) {
		ReadObjSize(file);
		std::ifstream fin(file.c_str());
		std::string order;
	
		int matID = -1;
		int vCnt = 0,  fCnt = 0;
		while (getline(fin, order, '\n')) {
			std::stringstream fin2(order);
			std::string var;
			if (!(fin2 >> var)) continue;
	
			
			if (var == "v") {
				vCnt++;
				v[vCnt].Input(fin2);
			}		
			if (var == "f") {
				Triangle* tri = tris[fCnt] = new Triangle;
				tri->SetParent(polyhedron);
			
				if (matID != -1)
					tri->SetMaterial(mat[matID]);
				else
					tri->SetMaterial(polyhedron->GetMaterial());
				std::string str;
				for (int i = 0; fin2 >> str; i++) {
					int bufferLen = 0, buffer[3];
					buffer[0] = buffer[1] = buffer[2] = -1;
					for (int s = 0, t = 0; t < (int)str.length(); t++)
						if (t + 1 >= (int)str.length() || str[t + 1] == '/') {
							buffer[bufferLen++] = atoi(str.substr(s, t - s + 1).c_str());
							s = t + 2;
						}
					int vertexID = i;
					if (i >= 3) {
						vertexID = 2;
						tri = tris[fCnt] = new Triangle;
						*tri = *tris[fCnt - 1];
						tri->GetVertex(1) = tri->GetVertex(2);
						tri->GetPos(1) = tri->GetPos(2);
						tri->GetTextureVertex(1) = tri->GetTextureVertex(2);
						tri->GetNormalVectorID(1) = tri->GetNormalVectorID(2);
					}
					if (buffer[0] > 0) {
						tri->GetVertex(vertexID) = buffer[0];
						VDouble3 vertexPos = v[buffer[0]];
						vertexPos = vertexPos.Rotate(VDouble3(1, 0, 0), polyhedron->GetAngles().GetCoord(0));
						vertexPos = vertexPos.Rotate(VDouble3(0, 1, 0), polyhedron->GetAngles().GetCoord(1));
						vertexPos = vertexPos.Rotate(VDouble3(0, 0, 1), polyhedron->GetAngles().GetCoord(2));
						vertexPos = polyhedron->GetO() + vertexPos * polyhedron->GetSize();
						tri->GetPos(vertexID) = vertexPos;
					}
					if (buffer[1] > 0) {
						tri->GetTextureVertex(vertexID) = buffer[1];
					}
					if (buffer[2] > 0) {
						tri->GetNormalVectorID(vertexID) = buffer[2];
					}
					if (i >= 2) {
						tri->PreTreatment();
						fCnt++;
					}
				}
			}
		}
		fin.close();		
		CalnVn();	
		TriangleNode* root = polyhedron->GetTree()->GetRoot();
		root->size = fCnt;
		root->tris = new Triangle*[root->size];
		for (int i = 0; i < root->size; i++) {
			root->tris[i] = tris[i];
			root->box.Update(tris[i]);
		}
		polyhedron->GetTree()->BuildTree();
		
		polyhedron->SetVertexN(vn);
		polyhedron->SetPixel(vt);
		delete[] v;
	
		if(fSize>0){
			//for ( int i = 0 ; i <fSize ; i++ ) delete  tris[i];
	 		delete[] tris;
	 		fSize=0;
		}
		if(matSize>0){
			//for ( int i = 0 ; i <matSize ; i++ ) delete  mat[i];
	 		delete[] mat;
	 		matSize=0;
		}
	}	
};



class Sphere : public NoLightObj {
	VDouble3 O , De , Dc;//O球心坐标,dedc球的坐标轴（z轴和与之垂直的辐角为0的轴），用于计算纹理
	double R;//球的半径

public:
	Sphere(VDouble3 O1 , VDouble3 De1 ,VDouble3 Dc1,double R1,ColorRGB color1 , ColorRGB absor1,
		double refl1 , double refr1,double diff1 , double spec1,
		double rindex1,double drefl1,std::string texture1 ):NoLightObj( color1 ,  absor1,
		 refl1 ,  refr1, diff1 ,  spec1,
		 rindex1, drefl1, texture1 	){ 
		 
		O=O1 , De=De1 , Dc=Dc1;
		R=R1;
	}
	~Sphere() {}

	ColliderRecord CollideOne( Ray ray ) {
		ColliderRecord colli_R;
		ray.dir = ray.dir.GetUnitVect();//?????????
		VDouble3 eo = O - ray.pos;//光源指向球心的向量
		double eA = eo .Dot( ray.dir );//光源到投影点的距离 eA= |eo|cosab；
		if (eA <= EPS) return colli_R;//1）如果光源在球体外部并且 eA < 0，那么光 线与球面不相交  
		
		double AB2 = R * R- eo.Module2() + eA * eA;//AB2投影点  到 光线与球面的交点的距离平方 	
		if (AB2 <= EPS) return colli_R;//2）AB<=0,即球心到光线的距离大于半径，那么光线与球面不相交 
		
			AB2 = sqrt( AB2 );
			double eB1 = eA - AB2  , eB2 = eA + AB2; 
			if ( eB2 < EPS ) return colli_R;
			 
			if ( eB1 > EPS ) {
				colli_R.dista = eB1;
				colli_R.frontage = true;
			} else {
				colli_R.dista = eB2;
				colli_R.frontage = false;
			} 
		colli_R.impact = true;
		colli_R.In = ray.dir;
		colli_R.priBs=this;
		colli_R.Co = ray.pos + ray.dir * colli_R.dista;
		colli_R.No = ( colli_R.Co - O ).GetUnitVect();
		if (!colli_R.frontage) colli_R.No = -colli_R.No;
		return colli_R;
	}
	ColorRGB GetTexture(VDouble3 C) {
		if (material->texture != NULL)return ColorRGB(1,1,1);
		VDouble3 I = ( C - O ).GetUnitVect();
		double a = acos( -I.Dot( De ) );
		double b = acos( std::min( std::max( I.Dot( Dc ) / sin( a ) , -1.0 ) , 1.0 ) );
		double u = a / PI , v = b / 2 / PI;
		if ( I.Dot(Dc.Cross(De)) < 0 ) v = 1 - v;
		return material->texture->GetSmoothColor( u , v );
	}
};

class Plane : public NoLightObj {
	VDouble3 N , Dx , Dy;//N平面法向量
	//dx,dy平面的坐标轴，用于计算纹理，纹理图片在场景中的长宽
	double R;//平面与原点距离

public://Plane(): NoLightObj() { }
	Plane(VDouble3 N1 , VDouble3 Dx1 ,VDouble3 Dy1,double R1,ColorRGB color1 , ColorRGB absor1,
		double refl1 , double refr1,double diff1 , double spec1,
		double rindex1,double drefl1,std::string texture1 ):NoLightObj( color1 ,  absor1,
		 refl1 ,  refr1, diff1 ,  spec1,
		 rindex1, drefl1, texture1 	){ 
		 
		N=N1 , Dx=Dx1 , Dy=Dy1;
		R=R1;
		N = N.GetUnitVect();
	}

	~Plane() {}

	ColliderRecord CollideOne( Ray ray ) {
		ColliderRecord colli_R;
		ray.dir = ray.dir.GetUnitVect();
		N = N.GetUnitVect();
		double nomDi = N.Dot( ray.dir );
		if ( fabs( nomDi ) <= EPS ) return colli_R;//Point1 =N*R---N*Point1=R
		double tt = (  R-ray.pos.Dot( N ) )/ nomDi;// ( N * R - ray.pos ).Dot( N ) / nomDi;
		if ( tt <= EPS ) return colli_R;//(Point1 - ray.Position) * normal / eo;
		
		colli_R.impact = true;
		colli_R.In = ray.dir;
		colli_R.priBs=this;
		colli_R.dista = tt;
		colli_R.frontage = ( nomDi < 0 );
		colli_R.Co = ray.pos + ray.dir * colli_R.dista;
		colli_R.No = ( colli_R.frontage ) ? N : -N;
		
		return colli_R;
	}
	ColorRGB GetTexture(VDouble3 C) {
		if (material->texture != NULL)return ColorRGB(1,1,1);
		double u = C.Dot( Dx ) / Dx.Module2() + 0.5;
		double v = C.Dot( Dy ) / Dy.Module2() + 0.5;
	
		return material->texture->GetSmoothColor( u , v );
	}
};

class AreaLight : public LightObject {
	VDouble3 O , Dx , Dy;//O光源中心的位置 	//,dx矩形光源的x半轴,dy矩形光源的y半轴,光源的坐标，其长度有意义（半轴长）
	VDouble3 N ;//法向量 
	double R;//原点到面光源平面的距离 
public:
	AreaLight() : LightObject() {}
	AreaLight(VDouble3 O1,VDouble3 Dx1,VDouble3 Dy1,ColorRGB rgb1)  {
		O=O1;
		Dx=Dx1;
		Dy=Dy1; //,dx矩形光源的x半轴,dy矩形光源的y半轴,光源的坐标，其长度有意义（半轴长）
 		color=rgb1;
		
		N = (Dx.Cross(Dy)).GetUnitVect();
	 	R=O.Dot( N ) ;
	}
	~AreaLight() {}	

	ColliderRecord CollideOne(  Ray ray ) {
		ColliderRecord ret;
		
		ret.priBs=this;	
		ray.dir = ray.dir.GetUnitVect();//VDouble3 N = (Dx.Cross(Dy)).GetUnitVect();
		double nomDi = N.Dot( ray.dir );
		if ( fabs( nomDi ) <= EPS ) return ret;
		double tt = (  R-ray.pos.Dot( N ) )/ nomDi;//double tt = ( N * R - ray.pos ).Dot( N ) / nomDi;
		if ( tt < EPS ) return ret;//这些事判断直线、平面相交 
	
		VDouble3 C = ( ray.pos + ray.dir * tt ) - O;//计算OC--下面事判断交点在矩形内 
		if ( fabs( Dx.Dot( C ) ) > Dx.Dot( Dx ) ) return ret;//OC在Dx的投影大于Dx则不相交 
		if ( fabs( Dy.Dot( C ) ) > Dy.Dot( Dy ) ) return ret;//OC在Dy的投影大于Dy则不相交 
	
		ret.impact = true;
		ret.dista = tt;
		return ret;
	}
	ColorRGB getIrra( ColliderRecord* colli_R , NoLightObj* NObj_head , int shade_quality , int* hash ) {
	//计算碰撞点的光照度乘上diff或spec的比率，考虑阴影 
		ColorRGB ret;
		Ray ray( colli_R->Co);
		for ( int i = -2 ; i < 2 ; i++ )//-2,-1,0,1,共4点 
			for ( int j = -2 ; j < 2 ; j++ )//-2,-1,0,1,共4点 
				for ( int k = 0 ; k < shade_quality ; k++ ) {//shade_quality= 4
					ray.dir = (O + Dx * ( ( ran() + i ) / 2 ) + Dy * ( ( ran() + j ) / 2 )) - colli_R->Co;//方向随机扰动 
					double dist = ray.dir.Module();
					bool shade = false;
					for(NoLightObj* nowobj = NObj_head ; nowobj != NULL ; nowobj =(NoLightObj*) nowobj->GetNext() ) {
						ColliderRecord colli_R = nowobj->CollideOne( ray);
						if ( colli_R.impact && colli_R.dista < dist ) {
							shade = true;
							break;
						}
					}	
					if ( shade == false )
	ret += color * (colli_R->priBs->GetMaterial()->BRDF(ray.dir, colli_R->No, -colli_R->In) /(dist*dist));
				}
		ret /= 16.0 * shade_quality;
		return ret;
	}
	
	ColorRGB CalnIrradiance( ColliderRecord* colli_R , VDouble3 V , double locadist2 ) {
		ColorRGB ret = color * colli_R->priBs->GetMaterial()->BRDF(V, colli_R->No, -colli_R->In)/locadist2;
		return ret;
	}
	Ray getPhoton(){
		Ray ret;
		ret.power = color / color.Average3();
		ret.pos = O + Dx * ( ran() * 2 - 1 ) + Dy * ( ran() * 2 - 1 );
		ret.dir.AssRandomVect();
		return ret;
	}
};


class Camera {
public:	
	VDouble3 O , N , Dx , Dy,Dx1,Dy1;//感光点位置O、摄像机的法向N、镜头长宽半轴Dx，Dy 
	int W , H;//照片的象素长（或宽）---  bmp图像的象素高
	double lens_W , lens_H;//镜头的长（或宽）与感光点到镜头距离的比值 
	double dofSample, aperture, focalLen;//景深效果的
	
	double shade_quality,drefl_quality;//计算软阴影	
	int max_photons,emit_photons;//PM时最大光子数,总发射光子数
	int sample_photons;//采样光子数---采样时所需光子数
	double sample_dist;//采样光子半径--采样时最大半径


	ColorRGB** data;//储存渲染出来的图片（每个象素的色彩）

	Camera() {	}
	Camera(VDouble3 O1,VDouble3 N1, 	int W1 , int H1,
		double lens_W1, double lens_H1, double dofSample1, double aperture1, double focalLen1
		,double shade_quality1,double drefl_quality1,  int max_photons1,int emit_photons1,
		int sample_photons1, double sample_dist1 ) {
		O=O1; 	N=N1;	
		W=W1;	H=H1;
		lens_W=lens_W1,lens_H=lens_H1;
		dofSample=dofSample1,aperture=aperture1,focalLen=focalLen1;	//dofSample=12;	aperture=0.15;	focalLen=8;//景深 
		
		shade_quality=shade_quality1, drefl_quality=drefl_quality1;
		max_photons=max_photons1,emit_photons=emit_photons1;
		sample_photons=sample_photons1,sample_dist=sample_dist1;
		
		N = N.GetUnitVect();
		Dx = N.GetVerticalVect();
		Dy = Dx.Cross(N) .GetUnitVect();
		Dx1=Dx*aperture;
		Dy1=Dy*aperture;
		Dx = Dx * lens_W ;
		Dy = Dy * lens_H ;
		
		data = new ColorRGB*[H];
		for ( int i = 0 ; i < H ; i++ )data[i] = new ColorRGB[W];
	}
	
	~Camera(){
		if ( data == NULL ) {
			for ( int i = 0 ; i < H ; i++ )  delete[] data[i];
			delete[] data;
			data=NULL;
		}
	}
	
	VDouble3 GetO() { return O; }
	VDouble3 GetN() { return N; }
	int GetDofSample() { return dofSample; }
	double GetAperture() { return aperture; }
	double GetFocalLen() { return focalLen; }
	int GetW() { return W; }
	int GetH() { return H; }
	ColorRGB GetColor255(int i, int j) { return data[i][j]; }
	void SetColor( int i , int j , ColorRGB color ) { data[i][j] = color; }
	double GetShadeQuality() { return shade_quality; }
	double GetDreflQuality() { return drefl_quality; }


	
	int GetMaxPhotons() { return max_photons; }
	int GetEmitPhotons() { return emit_photons; }
	int GetSamplePhotons() { return sample_photons; }
	double GetSampleDist() { return sample_dist; }


	VDouble3 Emit( double i , double j ) {//得到象素(i,j)对应的射出光线
		return N  + Dx * ( j /(W-1) - 0.5 )+ Dy * (  i / (H-1) - 0.5 );
		//return N + Dx* lens_W * (j / (W - 1) - 0.5)+Dy * lens_H * (i / (H - 1) - 0.5) + ;
	}

	void DofEmit(double i, double j, VDouble3* dof_O, VDouble3* dof_V) {
		//需要渲染景深时象素(i,j)对应的一个随机的射出光线（改变dof_O和dof_V）
		*dof_O = O + Dx1* (ran() * 2 - 1) + Dy1* (ran() * 2 - 1);//*dof_O = O + Dx * aperture * x + Dy * aperture * y;//*dof_V =  O + Emit(i, j) * focalLen - *dof_O;	
		*dof_V = O + Emit(i, j) * focalLen - *dof_O;//.GetUnitVect()
	}
	
	void Output( Bitmap* bmp ) {//把得到的色光数据导入到bmp文件中
		bmp->Initialize( H , W );
		for ( int i = 0 ; i < H ; i++ )
			for ( int j = 0 ; j < W ; j++ )
				bmp->SetColor( i , j , data[i][j] );
	}
};

class Scene {
public:   
	NoLightObj* NObj_head;//物品集的链表头
	LightObject* LObj_head;//光源集的链表头
	Camera* camera;//摄像机的指针
	ColorRGB background_color;//背景颜色
	


	Scene() {
		NObj_head = NULL;
		LObj_head = NULL;
		background_color = ColorRGB();
	}
	~Scene(){
		while ( NObj_head != NULL ) {
			NoLightObj* next_head = (NoLightObj*)NObj_head->GetNext();
			delete NObj_head;
			NObj_head = next_head;
		}
	
		while ( LObj_head != NULL ) {
			LightObject* next_head = (LightObject*)LObj_head->GetNext();
			delete LObj_head;
			LObj_head = next_head;
		}
	
		delete camera;
	}
	
	NoLightObj* GetPrimitiveHead() { return NObj_head; }
	LightObject* GetLightHead() { return LObj_head; }
	Camera* GetCamera() { return camera; }
	ColorRGB GetBackgroundColor() { return background_color; }



	void AddPrimitive ( NoLightObj* pri ) {
		if ( pri != NULL ) {
			pri->SetNext( NObj_head );
			NObj_head = pri;
		}
	}
	
	void AddLight( LightObject* light ) {
		if ( light != NULL ) {
			light->SetNext( LObj_head );
			LObj_head = light;
		}
	}
	ColliderRecord* FindNearestCollide( Ray ray ) {
		ColliderRecord* ret = NULL;	//计算光线(ray_O,ray_V)碰到的第一个物体，以碰撞器类Collider包装返回
		for(NoLightObj* nowobj = NObj_head ; nowobj != NULL ; nowobj =(NoLightObj*) nowobj->GetNext() ) {
			ColliderRecord colli_R = nowobj->CollideOne( ray);
			if ( colli_R.impact && ( ret == NULL || colli_R.dista < ret->dista ) ) {
				if (ret == NULL) ret = new ColliderRecord;
				*ret = colli_R;
			}
		}
		return ret;
	}
	ColliderRecord* FindNearestLight( Ray ray ) {
		ColliderRecord* ret = NULL;	//计算光线(ray_O,ray_V)碰到的第一个光源，以光源碰撞器类ColliderRecord包装返回
		for(LightObject*nowobj = LObj_head ; nowobj != NULL ; nowobj =(LightObject*) nowobj->GetNext() ) {
			ColliderRecord lightCollider = nowobj->CollideOne( ray);
			if (lightCollider.impact && (ret == NULL || lightCollider.dista < ret->dista)) {
				if (ret == NULL) ret = new ColliderRecord;
				*ret = lightCollider;
			}
		}	
		return ret;
	}
	void CreateScene( ) {//从文件中读入场景的数据
	
		NoLightObj* new_primitive = NULL;
		LightObject* new_light = NULL;
		
		background_color=ColorRGB(0.1, 0.1, 0.1);
		
		camera=new Camera(VDouble3(0, -8 ,2),//视点坐标 
				VDouble3(0, 1, 0),//法向量 
				16,12,//照片的象素长w（高h）---  bmp图像的象素高 
				0.75,1.0,//镜头的长w（高h）与感光点到镜头距离的比值		
				32.0,0.15,8,//景深 //dofSample=128;	aperture=0.15;	focalLen=8
				//0.0,  0.0,  0.0,//景深 //dofSample=128;	aperture=0.15;	focalLen=8
				4.0,15.0,//计算软阴影  shade_quality1,  drefl_quality1		
				10000000,	//PM时最大光子数
				15000000,	//总发射光子数	
				500	,0.1	//PM时碰撞点的采样光子数,采样时最大半径
					);	
		new_light = new AreaLight(VDouble3(0, 0, 6),//矩形光源的中心点
				VDouble3(0.5, 0, 0),//矩形光源的x半轴
				VDouble3(0, 0.5, 0),//矩形光源的y半轴
				ColorRGB(50, 50, 50)//光源颜色 
				);
		AddLight( new_light );
		
		
		new_primitive = new Plane(VDouble3(0, 0, 1),//法线 
				VDouble3(20, 0, 0),//矩形的x半轴
				VDouble3(0, 20, 0),0.0,//矩形的y半轴 //原点到平面的距离 
				ColorRGB(1,1,1),ColorRGB(),// 表面颜色， 吸收色
				0.3, 0.0, 0.5, 0.0,//反射系数，折射系数，漫反射系数，镜面反射 
				0.0,0.15,"floor.bmp"  //两种介质折射率比值,偏差半径,纹理bmp文件名 
				);
		AddPrimitive( new_primitive );//地板 
		new_primitive->PreTreatment();
		
		new_primitive = new Plane(VDouble3(0, 1, 0),//法线 
				VDouble3(20, 0, 0),//矩形的x半轴
				VDouble3(0, 0, 20),6.0,//矩形的y半轴 //原点到平面的距离 
				ColorRGB(1,1,1),ColorRGB(),//表面颜色， 吸收色
				0.0, 0.0, 0.5, 0.0,//反射系数，折射系数，漫反射系数，镜面反射 
				0.0,0.15,"marble.bmp"  //两种介质折射率比值,偏差半径,纹理bmp文件名 
				);
		AddPrimitive( new_primitive );//前面的墙 
		new_primitive->PreTreatment();
		
		
		new_primitive = new Sphere(VDouble3(-2, 0,0.99 ),//球心 
				VDouble3(0, 0, 1),//矩形的x半轴
				VDouble3(0, 1, 0),0.99,//矩形的y半轴 //原点到平面的距离 
				ColorRGB(1,0.5,0.5),ColorRGB(0,0,0),//表面颜色， 吸收色
				0.3, 0.0, 0.5, 0.6,//反射系数，折射系数，漫反射系数，镜面反射 
				0.0,0.25,"marble2.bmp"  //两种介质折射率比值,偏差半径,纹理bmp文件名 
				);
		AddPrimitive( new_primitive );//个大理石球
		new_primitive->PreTreatment();
		
		
		new_primitive = new Polyhedron( VDouble3(0, 3,2.13),//球心 
				VDouble3(6, 6, 6),//size,龙的最大外围尺寸 
				VDouble3(0, 0, 90),0.0,//物体的方向 
				ColorRGB(0.5,1,0.5),ColorRGB(0,0,0),//表面颜色， 吸收色
				0.0, 0.0, 0.5, 0.5,//反射系数，折射系数，漫反射系数，镜面反射 
				0.0,0.0,"dragon.obj"  //两种介质折射率比值,偏差半径,纹理bmp文件名 
				);
		AddPrimitive( new_primitive );//dragon.obj
		new_primitive->PreTreatment();
		
		
		
		new_primitive = new Sphere(VDouble3(0.5, -0.5, 1.01),//球心 
				VDouble3(1, 0, 0),//矩形的x半轴
				VDouble3(0, 0, 1),1,//矩形的y半轴 //原点到平面的距离 
				ColorRGB(1,1,1),ColorRGB(0,0,1),//表面颜色， 吸收色
				0.5, 1.0, 0.1, 0.3,//反射系数，折射系数，漫反射系数，镜面反射 
				0,0 ,""  //两种介质折射率比值,偏差半径,纹理bmp文件名 
				);
		AddPrimitive( new_primitive );//一个透明的0.4玻璃球
		new_primitive->PreTreatment();
	}


};





	
	
	
	
	
	
Polyhedron::Polyhedron(VDouble3 O1 , VDouble3 size1 ,VDouble3 angles1,double R1,ColorRGB color1 , ColorRGB absor1,
		double refl1 , double refr1,double diff1 , double spec1,
		double rindex1,double drefl1,std::string texture1 ):NoLightObj( color1 ,  absor1,
		 refl1 ,  refr1, diff1 ,  spec1,
		 rindex1, drefl1, "" 	){
		 
		O=O1;
		size=size1;	
		angles=angles1;
		meshFile=texture1;
	
		angles *= PI / 180;		
		tree = new TriangleTree();
	}
	
	Polyhedron::~Polyhedron() {
		delete tree;
		if (vertexN != NULL)
			delete[] vertexN;
		if (pixel != NULL)
			delete[] pixel;
	}
	

	void Polyhedron::PreTreatment(){
		ObjReader* objReader = new ObjReader();
		objReader->SetPolyhedron(this);
		objReader->ReadObj(meshFile);
		delete objReader;
	}
	ColliderRecord Polyhedron::CollideOne(Ray ray) {
		return tree->CollideOne( ray);
	}
	
	
void PhotonTracing( Ray , int dep , bool refracted ); 
ColorRGB RayTracing(Ray ray, int dep, bool refracted, int* hash, int rc, ColorRGB weight);	
 
	const int MAX_PHOTONTRACING_DEP = 20;
	double start_time ;
	
	Scene* scene;//待渲染的场景 
 	int hh,ww;
	LightObject* light;
	double totalPower;
	
	const int MAX_DREFL_DEP = 2;
	const int MAX_RAYTRACING_DEP = 20;
	const int HASH_FAC = 7;
	const int HASH_MOD = 10000007;	

	Camera* camera;//场景中的摄像头
	int  H,W;
	int** hash_d;
	
	bool completeThread[10];
	std::mutex* completeRow;
	bool completeThread2[10];
	std::mutex* completeRow2;

	int max_photons , stored_photons , emit_photons;
	Ray* photons;
	VDouble3 box_min , box_max;//空间范围的两个角,包含所有的光子 
	
	
void SplitA( Ray* porg , int st , int en , int med , int axis ) {
	int l = st , r = en;
	while ( l < r ) {
		double key = porg[r].pos.GetCoord( axis );//以末尾元素为候选点 
		int i = l - 1 , j = r;
		for ( ; ; ) {
			while ( porg[++i].pos.GetCoord( axis ) < key );//从左边向中间扫描，得到较大的点 
			while ( porg[--j].pos.GetCoord( axis ) > key && j > l );//从右边向中间扫描 ，得到较小的点 
			if ( i >= j ) break;
			std::swap( porg[i] , porg[j] );
		}
		std::swap( porg[i] , porg[r] );
		if ( i >= med ) r = i - 1;
		if ( i <= med ) l = i + 1;
	}
}
void BalanceSegment( Ray* porg , int index , int st , int en ) {
	if ( st == en ) {//1）终止条件 
		photons[index] = porg[st];//创建叶子
		return;
	}
	int med = 1;//2）得到 光子个数的中间值 
	while ( 4 * med <= en - st + 1 ) med <<= 1;
	if ( 3 * med <= en - st + 1 ) med = med * 2 + st - 1;
		else med = en + 1 - med;	
	int axis = 2;//3）得到切分轴(或切分平面)，0=x，1=y，2=Z ---包围所有光子的两个角点（最大、最小），x、y、z三个轴向上差距最大的那个轴 
	if ( box_max.x - box_min.x > box_max.y - box_min.y && box_max.x - box_min.x > box_max.z - box_min.z ) axis = 0; 
	else if ( box_max.y - box_min.y > box_max.z - box_min.z ) axis = 1;
	SplitA( porg , st , en , med , axis );//得到切分中间值 
	photons[index] = porg[med];  
	photons[index].plane = axis;// 将数据集分为两个子集
	if ( st < med ) {
		double tmp = box_max.GetCoord( axis );
		box_max.GetCoord( axis ) = photons[index].pos.GetCoord( axis );
		BalanceSegment( porg , index * 2 , st , med - 1 ); //子递归构造左子树  集递归调用buildKdTree函数  
		box_max.GetCoord( axis ) = tmp;
	}
	if ( med < en ) {
		double tmp = box_min.GetCoord( axis );
		box_min.GetCoord( axis ) = photons[index].pos.GetCoord( axis );
		BalanceSegment( porg , index * 2 + 1 , med + 1 , en ); //递归构造右子树  子集递归调用buildKdTree函数  
		box_min.GetCoord( axis ) = tmp;
	}
}

void LocatePhotons( Near_pho* np , int index ) {
//	Kd树的查询：判断是否叶子节点，是则将叶子节点的温度加入累加和、累加个数增1；判断左子树完全在所查询区域内，则该节点的温度加入累加和、累加个数。判断左子树与所查询区域相交且不空，则递归查询左子树；判断右子树完全在所查询区域内，则该节点的温度加入累加和、累加个数。判断右子树与所查询区域相交且不空，则递归查询右子树
	if ( index > stored_photons ) return;
	Ray *photon = &photons[index];

	if ( index * 2 <= stored_photons ) {
		double dist = np->pos.GetCoord( photon->plane ) - photon->pos.GetCoord( photon->plane );
		if ( dist < 0 ) {
			LocatePhotons( np , index * 2 );
			if ( dist * dist < np->dist2[0] ) LocatePhotons( np , index * 2 + 1 );
		} else {
			LocatePhotons( np , index * 2 + 1 );
			if ( dist * dist < np->dist2[0] ) LocatePhotons( np , index * 2 );
		}
	}

	double dist2 = photon->pos.Distance2( np->pos );
	if ( dist2 > np->dist2[0] ) return;

	if ( np->found < np->max_photons ) {
		np->found++;
		np->dist2[np->found] = dist2;
		np->photons[np->found] = photon;
	} else {
		if ( np->got_heap == false ) {
			for ( int i = np->found >> 1 ; i >= 1 ; i-- ) {
				int par = i;
				Ray* tmp_photon = np->photons[i];
				double tmp_dist2 = np->dist2[i];
				while ( ( par << 1 ) <= np->found ) {
					int j = par << 1;
					if ( j + 1 <= np->found && np->dist2[j] < np->dist2[j + 1] ) j++;
					if ( tmp_dist2 >= np->dist2[j] ) break;
					
					np->photons[par] = np->photons[j];
					np->dist2[par] = np->dist2[j];
					par = j;
				}
				np->photons[par] = tmp_photon;
				np->dist2[par] = tmp_dist2;
			}
			np->got_heap = true;
		}

		int par = 1;
		while ( ( par << 1 ) <= np->found ) {
			int j = par << 1;
			if ( j + 1 <= np->found && np->dist2[j] < np->dist2[j + 1] ) j++;
			if ( dist2 > np->dist2[j] ) break;

			np->photons[par] = np->photons[j];
			np->dist2[par] = np->dist2[j];
			par = j;
		}
		np->photons[par] = photon;
		np->dist2[par] = dist2;

		np->dist2[0] = np->dist2[1];
	}
}

ColorRGB getIrra( ColliderRecord* colli_R , double max_dist , int n ) {
	ColorRGB ret;

	Near_pho np;
	np.pos = colli_R->Co;   
	np.max_photons = n;
	np.dist2 = new double[n + 1];
	np.photons = new Ray*[n + 1];
	np.dist2[0] = max_dist * max_dist;

	LocatePhotons( &np , 1 );
	if ( np.found <= 8 ) return ret;

	for ( int i = 1 ; i <= np.found ; i++ ) {
		VDouble3 dir = np.photons[i]->dir;
		if ( colli_R->No.Dot( dir ) < 0 )
			ret += np.photons[i]->power * colli_R->priBs->GetMaterial()->BRDF(-dir, colli_R->No, -colli_R->In);
	}
	
	ret = ret * (4 / (emit_photons * np.dist2[0]));
	return ret;
}



void PhotonTracing(Ray photon, int dep, bool refracted) {
	if (dep > MAX_PHOTONTRACING_DEP) return;
	ColliderRecord* colli_R = scene->FindNearestCollide(Ray(photon.pos, photon.dir));
	if (colli_R == NULL) return;//没有碰撞，则递归停止 
	
	NoLightObj* pri = (NoLightObj*)colli_R->priBs;
	photon.pos = colli_R->Co;
	if (pri->GetMaterial()->diff > EPS && dep > 1) {
		if (stored_photons < max_photons) {
			photons[++stored_photons] = photon;
			box_min = VDouble3(std::min(box_min.x, photon.pos.x), std::min(box_min.y, photon.pos.y), std::min(box_min.z, photon.pos.z));
			box_max = VDouble3(std::max(box_max.x, photon.pos.x), std::max(box_max.y, photon.pos.y), std::max(box_max.z, photon.pos.z));
		}
	}

	double prob = 1;
	bool ret = true;
	MaterialQua* material = pri->GetMaterial();//物体材质 
	ColorRGB color = material->color;//物体颜色 
	if (material->texture != NULL) color = color * pri->GetTexture(colli_R->Co);//物体纹理 
	double eta = material->diff * color.Average3();//漫反射乘以颜色平均 
	if (eta <= ran() * (prob)) {//俄罗斯轮盘赌 
		prob -= eta;//prob变小了 
		ret = false;
	}else {
		photon.dir = colli_R->No.Diffuse();//碰撞点法线的漫反射光线 作为新光子方向 
		photon.power = photon.power * color / color.Average3();
		PhotonTracing(photon, dep + 1, refracted);//按照新方向、开始跟踪 
	}
	if (ret == false) {
		double eta = material->refl * color.Average3();//反射乘以颜色平均 
		if (eta <= ran() * (prob)) {//俄罗斯轮盘赌 
			prob -= material->refl;//prob变更小了 
			ret = false;
		}else {
			photon.dir = photon.dir.Reflect(colli_R->No);//碰撞点法线的反射光线 作为新光子方向 
			if (material->drefl > EPS) {//有偏差半径，就计算软阴影 
				VDouble3 Dx = photon.dir.GetVerticalVect();
				VDouble3 Dy = photon.dir.Cross(Dx).GetUnitVect();
				double dre = material->drefl*material->drefl;
				photon.dir += Dx * ((ran() * 2 - 1)*dre) + Dy * ((ran() * 2 - 1)*dre);//新光子方向 再增加一个偏差半径内的扰动 
			}
			photon.power = photon.power * color / color.Average3();
			PhotonTracing(photon, dep + 1, refracted);//按照新方向、开始跟踪  
			ret = true;
		}
		if (ret == false) {
			double eta = material->refr;//直接用折射系数做为光量 
			if (refracted) {//true表示在玻璃介质内部，向外（空气种界面） 
				ColorRGB trans = (material->absor * -colli_R->dista).Exp();//e的r次方，因r是负数， 
				eta *= trans.Average3();
				photon.power = photon.power * trans / trans.Average3();
			}
			if (eta <= ran() * (prob)) {//俄罗斯轮盘赌
				prob -= material->refr;//prob变更更小了 （第3次减少） 
			}else {
				double zhe = material->rindex;
				if (!refracted) zhe = 1 / zhe;
				bool nextRefracted = refracted;
				photon.dir = photon.dir.Refract(colli_R->No, zhe, &nextRefracted);
				PhotonTracing(photon, dep + 1, nextRefracted);
			}
		}
	}
	delete colli_R;

}

ColorRGB RayTracing(Ray ray , int dep , bool refracted , int* hash, int rc, ColorRGB weight) {
	if ( dep > MAX_RAYTRACING_DEP ) return ColorRGB();
	if ( hash != NULL ) *hash = ( *hash * HASH_FAC ) % HASH_MOD;
	ColorRGB ret;
	ColliderRecord* colli_R = scene->FindNearestCollide(ray );
	ColliderRecord* lightCollider = scene->FindNearestLight(ray );
	if (lightCollider != NULL) {
		LightObject* nearest_light = (LightObject*)lightCollider->priBs;
		if (colli_R == NULL || colli_R->dista < colli_R->dista) {
			if ( hash != NULL ) *hash = ( *hash + nearest_light->getHash_d() ) % HASH_MOD;
			ret += nearest_light->GetColor255() / nearest_light->GetColor255().RGBMax();
		}
		delete lightCollider;
	}
	if ( colli_R != NULL ) {
		NoLightObj* pri = (NoLightObj*)colli_R->priBs;//pri是最近的物体 
		if ( hash != NULL ) *hash = ( *hash + pri->getHash_d() ) % HASH_MOD;				
		if ( pri->GetMaterial()->diff > EPS ) {//漫反射  			
			ColorRGB color = pri->GetMaterial()->color;//材质颜色 
			if ( pri->GetMaterial()->texture != NULL ) color = color * pri->GetTexture(colli_R->Co);//材质纹理颜色 
			ret += color * scene->GetBackgroundColor() * pri->GetMaterial()->diff;//材质漫反射 
			for ( LightObject* light = scene->GetLightHead() ; light != NULL ; light =(LightObject*) light->GetNext() )//直接光--软阴影 
				ret += color * light->getIrra( colli_R , scene->GetPrimitiveHead() , scene->GetCamera()->GetShadeQuality() , hash );
			ret += color * getIrra( colli_R , camera->GetSampleDist() , camera->GetSamplePhotons() );//光子颜色 
		}
		if ( pri->GetMaterial()->refl > EPS ){
			VDouble3 dir = ray.dir.Reflect( colli_R->No );
			if ( pri->GetMaterial()->drefl < EPS || dep > MAX_DREFL_DEP ) {
				ColorRGB alpha = pri->GetMaterial()->color * pri->GetMaterial()->refl;
				ret += RayTracing(Ray(colli_R->Co,dir) , dep + 1 , refracted , hash, rc, weight * alpha) * alpha;
			}else{
				int totalSample = camera->GetDreflQuality();//阴影 
				ColorRGB ret1, alpha = pri->GetMaterial()->color * pri->GetMaterial()->refl / totalSample;
				VDouble3 Dx = dir.GetVerticalVect();	
				VDouble3 Dy = dir.Cross(Dx).GetUnitVect();
				double dre=pri->GetMaterial()->drefl*pri->GetMaterial()->drefl;
				for ( int k = 0 ; k < totalSample ; k++ ) {
					double x = ran() * 2 - 1,y = ran() * 2 - 1;			
					ret1 += RayTracing(Ray(colli_R->Co,  dir+Dx*(x*dre)+Dy*(y*dre)   ), dep + MAX_DREFL_DEP , refracted , NULL, rc, weight * alpha);
				}
				ret +=   ret1* alpha;
			}		
		}
		if ( pri->GetMaterial()->refr > EPS ) {
			double zhe = pri->GetMaterial()->rindex;
			if ( !refracted ) zhe = 1 / zhe;//false,表示在空气中，或者是在物体外面 
			bool nextRefracted = refracted;
			VDouble3 dir = ray.dir.Refract( colli_R->No , zhe , &nextRefracted );
			ColorRGB alpha = ColorRGB(1, 1, 1) * pri->GetMaterial()->refr;
			if (refracted)
				alpha *= (pri->GetMaterial()->absor * -colli_R->dista).Exp();//3个分量，分别求指数，e的r次方,e的g次方,e的b次方。
			ret += RayTracing(Ray(colli_R->Co,dir), dep + 1 , nextRefracted , hash, rc, weight * alpha)* alpha;
		}
		delete colli_R;
	}
	if ( dep == 1 ) ret = ret.AmLimit();
	return ret;
}



	int dofSample ;//dofSample= 128	
	float do1;//=1.0/(dofSample+1);
	

void Emitting(int threadID ) {
	int i=0;
	for (i = 0; i < hh ; i++) {
		if (!completeRow[i].try_lock()) continue;
		
		printf("A%d/%d\n", i, hh);	
		for ( int j = 0 ; j < ww ; j++ ) {
			Ray photon = light->getPhoton() ;
			photon.power *= totalPower;
			PhotonTracing( photon , 1 , false );		
		}
	}
	completeThread[threadID] = true;
}
void Sampling(int threadID ) {
	for ( int i = 0 ; i < H ; i++ ) {
		if (!completeRow[i].try_lock()) continue;
		printf("B%d/%d\n", i, H); 
		for ( int j = 0 ; j < W ; j++ ) {//48
			//if (!completeRow2[i].try_lock()) continue;
			//if (i%10 == 0) printf("2B%d/%d\n", i, H); 
			if (camera->GetAperture() < EPS) {	//aperture光圈大小		
				hash_d[i][j] = 0;
				ColorRGB color = camera->GetColor255(i, j);
				color += RayTracing(Ray(camera->GetO(),camera->Emit( i , j ) ), 1, false, &hash_d[i][j], 
						i * W + j, ColorRGB(1, 1, 1));
				camera->SetColor(i, j, color.AmLimit());
			} else {
				VDouble3 dof_O, dof_V;
				ColorRGB color = camera->GetColor255(i, j);
				for (int k = 0; k < dofSample; k++) {				
					camera->DofEmit(i, j, &dof_O, &dof_V);//DofEmit里面含有focalLen焦距 
					color+=RayTracing(Ray(dof_O, dof_V),1,false,NULL,i*W+j,ColorRGB(do1, do1, do1) );
				}
				camera->SetColor(i, j, (color/(dofSample+1)).AmLimit());
			}
		}
	}
	completeThread[threadID] = true;
}
void reSampling(int threadID ) {
	for ( int i = 0 ; i < H ; i++ ) {
		if (!completeRow[i].try_lock()) continue;
		printf("C%d/%d\n", i, H);
		for ( int j = 0 ; j < W ; j++ ) {//48
			if ( (i == 0 || hash_d[i][j] == hash_d[i - 1][j])&&(i == H - 1 || hash_d[i][j] == hash_d[i + 1][j])
			&&(j == 0 || hash_d[i][j] == hash_d[i][j - 1])&&(j == W - 1 || hash_d[i][j] == hash_d[i][j + 1]))continue;
			ColorRGB color = camera->GetColor255(i, j);
			for ( int r = -1 ; r <= 1 ; r++ )
				for ( int c = -1 ; c <= 1 ; c++ ) {
					if (((r + c) & 1) == 0) continue;
					color += RayTracing(Ray(camera->GetO(),  camera->Emit( i + ( double ) r / 3 , j + ( double ) c / 3 )
						) , 1, false, NULL, i * W + j, ColorRGB(0.2,0.2,0.2) );
				}
			camera->SetColor( i , j ,( color / 5).AmLimit() );
		}
	}
	completeThread[threadID] = true;
}
void MultiThreadSampling(int fg,int H1,void*pp ) {
	for (int i = 0; i < PM_MAX_THREADS; i++) completeThread[i] = false;
	completeRow = new std::mutex[H1];
	int i=0;
	for (i = 0; i < PM_MAX_THREADS; i++) {
		if(fg==0) {
			std::thread subThread(&Sampling,  i);
			if (i == PM_MAX_THREADS - 1) subThread.join();
			else subThread.detach();
		}else if(fg==1) {
			std::thread subThread(&reSampling,    i);
			if (i == PM_MAX_THREADS - 1) subThread.join();
			else subThread.detach();
		}else{//std::thread t1(&MyClass::func1,&myclass);  //std::thread subThread(&Photontracer::Emitting,  (Photontracer*)pp,  i);
			std::thread subThread(&Emitting,  i);
			if (i == PM_MAX_THREADS - 1) subThread.join();
			else subThread.detach();
		}
	} 
	for (bool end = false; !end; ) {
		end = true;
		for (int j = 0; j < PM_MAX_THREADS; j++){		
			if (!completeThread[j]) end = false;			
		}
	}
	delete[] completeRow;
}
int main() {	
	start_time = clock();	
	scene = new Scene;
	scene->CreateScene(  );
	camera = scene->GetCamera();
	dofSample = camera->GetDofSample();//dofSample= 128	
	do1=1.0/(dofSample+1);
	
	max_photons = scene->GetCamera()->GetMaxPhotons();
	stored_photons = 0;
	photons = new Ray[max_photons + 1];
	box_min = VDouble3( INF , INF , INF );
	box_max = VDouble3( -INF , -INF , -INF );
	emit_photons = scene->GetCamera()->GetEmitPhotons();		
	totalPower = 0;
	for (  light = scene->GetLightHead() ; light != NULL ; light =(LightObject*) light->GetNext() )
		totalPower += light->GetColor255().Average3();
	double photonPower = totalPower / scene->GetCamera()->GetEmitPhotons();

	for (  light = scene->GetLightHead() ; light != NULL ; light =(LightObject*) light->GetNext() ) {
		int lightPhotons = (int)(light->GetColor255().Average3() / photonPower);//lightPhotons = lightPhotons / PM_MAX_THREADS + ((lightPhotons % PM_MAX_THREADS > threadID) ? 1 : 0);
		hh=sqrt(lightPhotons) ;	
		ww=lightPhotons/hh;
		MultiThreadSampling( 2 ,hh,NULL);
		//for (int i = 0; i < hh ; i++) {
			//if (!completeRow[i].try_lock()) continue;				
			/*for ( int j = 0 ; j < lightPhotons ; j++ ) {
				if (j%50 == 0) printf("A%d/%d\n", j, lightPhotons);
				Ray photon = light->getPhoton() ;
				photon.power *= totalPower;
				PhotonTracing( photon , 1 , false );		
			}*/
		//}		
	}	
	Ray* porg = new Ray[stored_photons + 1];
	for ( int i = 0 ; i <= stored_photons ; ++i )porg[i] = photons[i];	
	BalanceSegment( porg , 1 , 1 , stored_photons );
	delete[] porg;
	
	H = camera->GetH(); W = camera->GetW();//rayQ= Ray(camera->GetO());
	hash_d = new int*[H];
	for (int i = 0; i < H; i++) hash_d[i] = new int[W];
	MultiThreadSampling(0, H, NULL);
	//#pragma omp parallel for schedule(dynamic, 1)
	/*for ( int i = 0 ; i < H ; i++ ) {		
		//if (!completeRow[i].try_lock()) continue;
		if (i%10 == 0) printf("B%d/%d\n", i, H);//std::cout << "Sampling:   " << i << "/" << H << std::endl;//
		for ( int j = 0 ; j < W ; j++ ) {
			if (camera->GetAperture() < EPS) {	//aperture光圈大小		
				hash_d[i][j] = 0;
				ColorRGB color = camera->GetColor255(i, j);
				color += RayTracing(Ray(camera->GetO(),camera->Emit( i , j ) ), 1, false, &hash_d[i][j], i * W + j, ColorRGB(1, 1, 1));
				camera->SetColor(i, j, color.AmLimit());
			} else {
				int dofSample = camera->GetDofSample();//dofSample= 128	
				ColorRGB color = camera->GetColor255(i, j);
				for (int k = 0; k < dofSample; k++) {
					VDouble3 dof_O, dof_V;
					camera->DofEmit(i, j, &dof_O, &dof_V);//DofEmit里面含有focalLen焦距 
					color += RayTracing(Ray(dof_O, dof_V), 1, false, NULL, i * W + j, ColorRGB(1, 1, 1) / dofSample) / dofSample;
				}
				camera->SetColor(i, j, color.AmLimit());
			}
		}
	}*/
	if (camera->GetAperture() < EPS) MultiThreadSampling(1, H, NULL);	
	for (int i = 0; i < H; i++) delete[] hash_d[i];
	delete[] hash_d;
	Bitmap* bmp = new Bitmap;
	camera->Output(bmp);
	bmp->Output("pic.bmp");
	delete bmp;
	delete[] photons;
	delete  scene;
	printf("Escaped time: %.5lf\n", (clock() - start_time) / CLOCKS_PER_SEC);
	return 0;		
}




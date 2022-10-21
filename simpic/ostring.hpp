/*  Our Own Local String Library, by PeterM */

#ifndef _OSTRING_H
#define _OSTRING_H

#include<stdio.h>

#include<vector>
//using namespace std;

#include "cstring"

class ostring
{
protected:
	unsigned short len;
	unsigned short size;
	char *d;

	inline  ostring *scopy(const ostring *source);
	inline  ostring *scopy(const char *source);

	inline  ostring *srealloc(int newsize);

public:

	//constructors and destructor
	ostring();
	~ostring() {len=size=0; delete[] d;};

	ostring(const ostring& copy);
	ostring(const char *copy);

	ostring(char c);
	ostring(int i)
	{
		len = 0; size = 0; d = 0;
		*this = i;
	}

	//conversion
	operator char * () { return d;};

	// assignment
	ostring& operator = (const ostring& two);
	ostring& operator = (const char *two);
	ostring& operator = (char c);
	ostring& operator = (int i);

	// comparison
	friend  int operator==(const ostring& a,const ostring& b);
	friend  int operator==(const ostring& a,const char *b);
	friend  int operator!= (const ostring& a,const ostring& b);
	friend  int operator!= (const ostring& a,const char *b);

	//  concatenation
	ostring& operator += (const ostring& two);
	ostring& operator += (const char *two);
	ostring& operator += (char c);
	ostring& operator += (int i);
	friend  ostring operator + (const ostring& a,const ostring& b);

	void append(const ostring &add);
	void append(const char add);
	void append(const char *add);

	// output
	const char *c_str() const {return d;}
	// mindgame: depreciate to use this fn.
	//           This fn is made since the argument of xgrafix function is
	//          not 'const char *' but 'char *'.
	char *c_str2() const {return d;}

	const char * operator ()() const {return d;}

	// length
	int length() { return len; } ;

	// searching
	inline int find(const char *lookfor,int startpos=0) const;
	inline int find(char lookfor,int startpos=0) const;
	inline int find(const ostring& str, int startpos=0) const;

	inline int findlast(char lookfor, int startpos=0) const;
	inline int charcount(char lookfor) const;

	bool contains(char c) const{ return (find(c)!=-1)? true:false;};
	bool contains(char *c) const { return (strstr(d,c))? true:false;};
	bool contains(ostring two) const { return (strstr(d,two.d))? true:false;};

	// getting part, from startpos to endpos, and returning the ostring.
	inline ostring substr(int startpos, int length);
	inline ostring substr(int startpos);

	// replacing
	inline ostring replace(char *from, char *to);

	//  bracket operator
	char& operator [](int index);

	// stripping parts of the string
	ostring & strip(char nuke_me);
	ostring & stripcomment(const char *comment_str);

        inline std::vector<ostring> partition(char divider);

	inline bool isnested(char lch, char rch);
	inline bool isfloat() const;
	inline bool isint();
	inline int ndigits(int i);

	inline void spaces(int number);
};

inline ostring * ostring::scopy(const char *source){
	if (d) delete[] d;
	d = new char[strlen(source)+1];
	if(source) strcpy(d,source);
	else d[0]='\0';
	len=size=strlen(source);
	return this;
}

inline ostring * ostring::scopy(const ostring *source){  
	if(this==source) return this;
	if (d) delete[] d;
	d = new char[(source->size)+1];
	if(source->d) strcpy(d,source->d);
	else d[0]='\0';
	size=source->size;
	len=source->len;
	return this;
}

inline ostring * ostring::srealloc(int newsize)
{
	char *tmp= new char[newsize+1];
	if(d)
	{
		strcpy(tmp, d);
		delete[] d;
	}	
	else tmp[0]='\0';
	d = tmp;
	size=newsize;
	len = strlen(d);
	return this;
}

inline ostring::ostring() { len = 0; size = 0; d = 0; }

inline ostring::ostring(const ostring &copy) 
{ 
	len = 0; size = 0; d = 0; scopy(&copy);
}

inline ostring::ostring(const char *copy) {len=0; size=0; d=0; scopy(copy);}
inline ostring::ostring(char c) { len=1;size=1;d=new char[2];d[0]=c;d[1]='\0';}

inline ostring & ostring::operator =(const ostring& two)
{
	return *(scopy(&two));
}

inline ostring & ostring::operator =(const char * two)
{
	return *(scopy(two));
}

inline ostring & ostring::operator =( char c)
{ len=1; size=2; if(d) delete[] d; d = new char[2];
d[0]=c;d[1]='\0';
return *this;
}

inline ostring & ostring::operator =(int i)
{
	len=ndigits(i); size=len+1; if(d) delete[] d; d = new char[size];
	sprintf(d, "%d", i);

	return *this;
}

inline void ostring::append(const ostring& add)
{ register int newlen;
if(size < (newlen = len + add.len)) srealloc(newlen);
len = newlen;
strcat(d,add.d);
}

inline void ostring::append(const char * add)
{ register int newlen;
if(size < (newlen=len + strlen(add))) srealloc(newlen);
len=newlen;
strcat(d,add);
}

inline void ostring::append(const char add)
{
	if(size < (len + 1)) srealloc(size+1);
	d[len]=add;
	d[len+1]='\0';
	len+=1;
}

inline ostring& ostring::operator +=  (const ostring& two) {
	append(two);
	return *this;
}

inline ostring& ostring::operator +=  (const char * two) {
	append(two);
	return *this;
}
inline ostring& ostring::operator +=  (char c) {
	append(c);
	return *this;
}

inline ostring& ostring::operator += (int i)
		{
	int lll;
	char *ccc;
	lll=ndigits(i);
	ccc = new char[lll+1];

	sprintf(ccc, "%d", i);
	append(ccc);

	delete [] ccc;

	return *this;
		}

inline ostring operator + (const ostring&a, const ostring& b)
{
	ostring r;
	r.d = new char[(r.len=a.len + b.len) +1];
	strcpy(r.d,a.d);
	strcat(r.d,b.d);
	r.size=r.len;
	return r;
}


inline int operator==(const ostring& a, const ostring& b) 
		{
	return (strcmp(a.d,b.d)==0)? 1:0;
		}

inline int operator==(const ostring& a, const char * b) 
		{
	return (strcmp(a.d,b)==0)? 1:0;
		}

inline int operator!=(const ostring& a, const ostring& b) 
		{
	return (strcmp(a.d,b.d)!=0)? 1:0;
		}

inline int operator!=(const ostring& a, const char * b) 
		{
	return (strcmp(a.d,b)!=0)? 1:0;
		}


/*  true substring search instead of search for an assemblage of characters */

inline int ostring::find(const char *lookfor, int startpos) const
{ register int i=startpos;
register int j;
register int jmax=strlen(lookfor);
register const char *c,*b,*e;
if(len==0) return -1;
for(c=&d[startpos];*c;c++,i++) {
	/* check for the entire lookfor string in d */
	for(b=lookfor,j=0,e=c;j<jmax&&((i+j)<len);j++,b++,e++)
		if(*b!=*e) break;
	if(j==jmax) return i;
}
return -1;
}

inline int ostring::find(char lookfor, int startpos) const
{ register char *c;
register int i;
for(c=&d[startpos],i=0;*c;c++,i++)
	if(lookfor==*c) return i;
return -1;
}

inline int ostring::find(const ostring &lookfor, int startpos) const
{ 
	return find(lookfor.d,startpos);
}

inline int ostring::findlast(char lookfor, int startpos) const
{
	int nv = 0;
	int retval = -1;
	while((nv = find(lookfor, nv)) > 0)
	{
		retval = nv;
	}
	return retval;
}

inline ostring ostring::replace(char *from, char *to)
{
	int startval =0;
	int nextval;
	ostring newstring("");
	while((nextval = find(from, startval)) != -1)
	{
		newstring += substr(startval, nextval-startval-1);
		newstring += to;
		startval = nextval + strlen(from);
	}
	if(newstring.length()) { *this = newstring; }

	return *this;
}

inline int ostring::charcount(char lookfor) const
{
	int o;
	int i, l;
	o = i = l = 0;

	while(i != -1)
	{
		i = find(lookfor, l);
		if(i != -1)
		{
			l += i + 1;
			o++;
		}
	}

	return o;
}

inline char & ostring::operator [] (int i)
{
	return d[i];
}

inline ostring& ostring::strip(char nuke_me)
{
	int i, k, l;

	for(i=0; i < (l = strlen(d)); i++)
	{
		if(d[i] == nuke_me)
		{
			for(k=i; k < l; k++)
			{
				d[k] = d[k+1];
			}
			i--;
		}
	}
	len= strlen(d);
	return (*this);
}

inline ostring& ostring::stripcomment(const char *comment_str)
{
	int i = find(comment_str, 0);
	if(i > -1)
	{
		d[i] = '\0';
	}
	len = strlen(d);
	return (*this);
}

inline ostring ostring::substr(int startpos, int length) {
	ostring r;
	register int i,j;

	if(len<1) return r;

	r.d = new char[(r.size=length)+1];
	for(i=startpos,j=0;j<length&&i<len;i++,j++) r.d[j]=d[i];
	r.d[j]='\0';
	r.len=strlen(r.d);
	return r;
}

inline ostring ostring::substr(int startpos)
{
	return substr(startpos, len);
}

inline std::vector<ostring> ostring::partition(char divider)
{
        std::vector<ostring> tvec;
	ostring tstr("");
	int pos;
	for(int i=0; i < length(); )
	{
		if((pos = find(divider, i)) >= 0)
		{
			tstr = substr(i, pos);						// find returns relative position of new divider;
														// changed "pos-i" to "pos"; JK 2018-01-27
			if (tstr != "") tvec.push_back(tstr);
			i += pos+1;									// i is absolute position in string, pos is relative; JK 2018-01-27
		}
		// mindgame: to avoid endless loop
		else
		{
			tstr = substr(i);
			if (tstr != "") tvec.push_back(tstr);
			break;
		}
	}

	return tvec;
}

inline bool ostring::isint()
{
	int i, l;
	i = 0;
	l = strlen(d);

	if(l == 0) { return false; }

	if(d[0] == '-')
	{
		i = 1;
	}
	if(i == l)
	{
		return false;
	}

	for(; i < l; i++)
	{
		if((d[i] < '0') || (d[i] > '9'))
		{
			break;
		}
	}
	if(i < l) { return false; }
	return true;
}

inline bool ostring::isfloat() const
{
	int i, j, k, l;
	i = 0;
	l = strlen(d);

	if(l == 0) { return false; }

	if(d[0] == '-')
	{
		i = 1;
	}

	if(i == l)
	{
		return false;
	}

	if(charcount('e') + charcount('E') > 1)
	{
		return false;
	}

	if(charcount('.') > 1)
	{
		return false;
	}

	j = find('.');
	k = find('e');

	if(k < 0)
	{
		k = find('E');
	}
	if(k >= 0)
	{
		if(j > k)
		{
			return false;
		}
	}

	if((j == i) || (k == i))
	{
		return false;
	}

	if(k == l-1)
	{
		return false;
	}

	for(; i < l; i++)
	{

		if((i == j) || (i == k)) { continue; }
		if(d[i] == '-')
		{
			if((i == k+1) && (i < l-1))
			{
				continue;
			}
		}
		if((d[i] < '0') || (d[i] > '9'))
		{
			break;
		}
	}
	if(i < l) { return false; }
	return true;
}

// isnested returns whether lch or rch fall between
// a matched set of lch and rch
// returns false if there is no matched set
inline bool ostring::isnested(char lch, char rch)
{
	int lpos, rpos;
	int lp2;

	lpos = find(lch);
	if(lpos < 0) { return false; }
	rpos = find(rch, lpos);
	if(rpos < 0) { return false; }

	lp2 = find(lch, lpos+1);
	if((lp2 < rpos) || (lp2 > 0))
	{
		return true;
	}
	return false;
}

inline int ostring::ndigits(int i)
{
	int nd=1;
	if(i == 0) { return 1; }
	if(i < 0) { i = -i; nd++; } // treat the minus sign as a digit

	while(i /= 10)
	{
		nd++;
	}
	return nd;
}

inline void ostring::spaces(int number)
{
	(*this) = "";
	for(int i=0; i < number; i++)
	{
		(*this) += ' ';
	}
}

#endif



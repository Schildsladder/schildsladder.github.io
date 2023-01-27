/*

Utils.js
========

Miscellaneous utilities

*/

Utils = {};

//	Randomly chosen unit vector
//	===========================

Utils.randomUnitVec =
function (posZ)
{
var theta=Math.acos(1.0-2.0*Math.random()), phi=(2.0*Math.PI)*Math.random();
var sinTheta=Math.sin(theta), cosTheta=Math.cos(theta);
return [
	sinTheta*Math.cos(phi),
	sinTheta*Math.sin(phi),
	posZ ? Math.abs(cosTheta): cosTheta
	];
}


//	Rotation by a given angle around a given axis
//	=============================================

//	The rotation matrix to rotate by an angle A
//	around the (unit) vector axis[0,1,2] is
//
//	axis[l]axis[k](1-cos A) + I cos A + (axis[i]sign{i,l,k}) sin A
//
//	(N.B. This looks counterclockwise, i.e. conventionally +ve for +ve A,
//	if you're looking in the direction axis[])
//
//	(1-(l+4-k)%3) gives sign{3-l-k,l,k}, where {3-l-k} = {0,1,2} \ {l,k}, l!=k

Utils.setRotation =
function (rotation, angle, axis, mag)
{
var sinA=Math.sin(angle);
var cosA=Math.cos(angle);
var magsq=mag*mag;
var noRotation = mag<=0.0 || angle==0.0;

for (var l=0;l<3;l++)
	for (var k=0;k<3;k++)
		rotation[l][k] = noRotation ? (l==k ? 1 : 0) : (
			(1.0-cosA)*axis[l]*axis[k]/magsq +
			((l==k) ? cosA : (1-(l+4-k)%3)*sinA*axis[3-l-k]/mag)
				);
}

//	Quicksorts
//	==========

//	Sort array a[][] by a[*][n]

Utils.qs =
function (n, a, i1, i2)
{
var low=a[i1][n], high=a[i2][n];
for (var i=i1;i<=i2;i++)
	{
	var ai=a[i][n];
	if (ai<low) low=ai;
	if (ai>high) high=ai;
	};
Utils.qs2(n,a,i1,i2,low,high);
};

//	Sort array a[][] by a[*][n]; we must specify lowest and highest values

Utils.qs2 =
function (n, a, i1, i2, low, high)
{
var t, med=(low+high)/2;
if (i2>i1 && med>low && med<high)
	{
	var h1=low, l1=high, v;
	var k1=i1, k2=i2;
	while (true)
		{
		while (k1<=i2 && (v=a[k1][n])<med)
			{
			if (v>h1) h1=v; k1++;
			};
		while (k2>=i1 && (v=a[k2][n])>=med)
			{
			if (v<l1) l1=v; k2--;
			};
		if (k1<k2) {t=a[k1]; a[k1]=a[k2]; a[k2]=t;} else break;
		};
	if (k1-i1>i2-k2) {Utils.qs2(n,a,i1,k1-1,low,h1); Utils.qs2(n,a,k2+1,i2,l1,high);}
	else {Utils.qs2(n,a,k2+1,i2,l1,high); Utils.qs2(n,a,i1,k1-1,low,h1);};
	};
};


//	Create multi-dimensional array

Utils.multiDim0 =
function (dims,ival,indx)
{
var result;
var d=dims[indx];
result = new Array(d);
if (indx==dims.length-1) for (var i=0;i<d;i++) result[i]=ival;
else for (var i=0;i<d;i++) result[i]=Utils.multiDim0(dims,ival,indx+1);
return result;
};


Utils.multiDim =
function (dims,ival)
{
return Utils.multiDim0(dims,ival,0);
};

//	Dot product

Utils.dot =
function(v,w)
{
var n=Math.min(v.length,w.length);
var res=0;
for (var i=0;i<n;i++) res+=v[i]*w[i];
return res;
};

Utils.dot2 =
function(v,w,n)
{
var res=0;
for (var i=0;i<n;i++) res+=v[i]*w[i];
return res;
};

//	Normalise a vector

Utils.normalise =
function (v)
{
var mag=Math.sqrt(Utils.dot(v,v));
if (mag>0) for (var i=0;i<v.length;i++) v[i]/=mag;
return mag;
};

//	Cross product

Utils.cross =
function (v, w)
{
return [v[1]*w[2]-v[2]*w[1], v[2]*w[0]-v[0]*w[2], v[0]*w[1]-v[1]*w[0]];
}

Utils.cross2 =
function (v, w, c)
{
c[0]=v[1]*w[2]-v[2]*w[1];
c[1]=v[2]*w[0]-v[0]*w[2];
c[2]=v[0]*w[1]-v[1]*w[0];
}

//	Transform a vector with a given matrix
//	======================================

//	First routine puts result back into original vector

Utils.temp=null;
Utils.tempD=-1;

Utils.transform1 =
function (matrix, vector)
{
var dim=matrix.length;
if (dim!=Utils.tempD)
	{
	Utils.temp=new Array(dim);
	Utils.tempD=dim;
	};
var tmp=Utils.temp;

for (var j=0;j<dim;j++)
	{
	tmp[j]=0.0;
	var mj=matrix[j];
	for (var k=0;k<dim;k++) tmp[j]+=mj[k]*vector[k];
	};
for (var j=0;j<dim;j++) vector[j]=tmp[j];
}

//	Second routine puts result in a different vector

Utils.transform2 =
function (matrix, src, dest)
{
var dim=matrix.length;
for (var j=0;j<dim;j++)
	{
	dest[j]=0.0;
	var mj=matrix[j];
	for (var k=0;k<dim;k++) dest[j]+=mj[k]*src[k];
	};
}

//	Test for infinite value

Utils.isInfinite =
function (f)
{
return !isFinite(f);
};

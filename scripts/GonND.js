/*

GonND.js
========

Class representing a planar polygon in n dimensions.

It is up to the calling program to provide vertices that are all COPLANAR, and
which do NOT have any edge-edge intersections.  A GonND itself can determine
the convexity/nonconvexity of the polygon.

GonND implements the interface "BSPobject", which defines a method by which an
object can split itself into parts for a Binary Space Partition (class "BSP").

*/


/*

Original Java class elements:

public int dim;						//	Dimension of space
public int nvert;					//	Number of vertices
public double vertices;			//	List of vertex (x,y,z,...) coordinates

private int convexity;				//	Convexity:  0=no, 1=yes, -1=don't know

public int surfaceID;				//	A user-defined surface ID
									//	Inherited by any GonNDs created by splitting.

public int flags;					//	User-defined flags
									//	Inherited by any GonNDs created by splitting.
								
public int edgeIDs;				//	If not null, an array of user-defined edge IDs.
									//	Inherited by any GonNDs created by splitting.
									
*/
									

//	Constructor
//	===========

//	(1) Constructor to form a polygon by simply copying a pointer to a complete array
//	of vertices.
//	(2) Constructor to form a polygon by simply copying a pointer to a complete array
//	of vertices; convexity is known.
//	(3) Constructor to form a polygon by taking copies of the individual vertex pointers.
//	(4) Constructor to form a polygon by taking copies of the data; convexity is known.

function GonND(d, v)
{
this.dim=d;
var narg=arguments.length;

if (narg==3)
	{
	//	(1) GonND(dim, vertices, nvert)
	this.vertices=v;
	this.nvert=arguments[2];
	this.convexity=this.nvert<=3 ? 1 : -1;
	}
else if (narg==4)
	{
	var boolArg = arguments[3]===true || arguments[3]===false;
	if (boolArg)
		{
		//	(2) GonND(dim, vertices, nvert, isConvex)
		this.vertices=v;
		this.nvert=arguments[2];
		this.convexity=arguments[3] ? 1 : 0;
		}
	else
		{
		//	(3) GonND(dim, varray, start, nvert)
		var start=arguments[2];
		this.nvert=arguments[3];
		this.vertices=new Array(this.nvert);
		for (var i=0;i<this.nvert;i++) this.vertices[i]=v[start+i];
		this.convexity=this.nvert<=3 ? 1 : -1;
		};
	}
else if (narg==5)
	{
	//	(4) GonND(dim, varray, start, nvert, isConvex)
	var start=arguments[2];
	this.nvert=arguments[3];
	this.vertices=new Array(this.nvert);
	for (var i=0;i<this.nvert;i++) this.vertices[i]=v[start+i];
	this.convexity=arguments[4] ? 1 : 0;
	}
else
	{
	throw "Bad constructor for GonND"
	};
	
this.surfaceID=0;
this.flags=0;
this.edgeIDs=null;
}

//	Flags and ID field for edgeIDs

GonND.EDGE_ORIGINAL=(1<<30);		//	Original unsplit edge
GonND.EDGE_INTERNAL=(1<<29);		//	Fake edge for regions
									//	with multiple boundaries
									//	(e.g. annular regions)
GonND.EDGE_POSITIVE=(1<<28);		//	Piece on +ve side of split
GonND.EDGE_NEGATIVE=(1<<27);		//	Piece on -ve side of split
GonND.EDGE_FROM_SPLIT=GonND.EDGE_POSITIVE|GonND.EDGE_NEGATIVE;
GonND.EDGE_IN_CONTENT=(1<<26);		//	This edge segment inside
									//	contents of splitting plane
GonND.EDGE_CUT_AT_START=(1<<25);	//	Set (along with EDGE_ORIGINAL)
GonND.EDGE_CUT_AT_END=(1<<24);		//	when original edge is not whole
GonND.EDGE_CUT=GonND.EDGE_CUT_AT_START|GonND.EDGE_CUT_AT_END;
													
GonND.EDGE_ID=(1<<16)-1;			//	Bit mask for edge ID number


//	newGon()
//	========

//	newGon() returns a new GonND by copying vertex pointers;
//	subclasses override this, to return a GonND with the full specificity of
//	the subclass.

GonND.prototype.newGon =
function (varray, start, nv)
{
return new GonND(this.dim, varray, start, nv);
}
	
//	translate()
//	===========

GonND.prototype.translate =
function (delta)
{
for (var i=0;i<this.nvert;i++)
	{
	var vi=this.vertices[i];
	for (var j=0;j<this.dim;j++) vi[j]+=delta[j];
	};
};

//	com()
//	=====

//	Centre of mass

GonND.prototype.com =
function()
{
var c=new Array(this.dim);
for (var j=0;j<this.dim;j++) c[j]=0;
for (var i=0;i<this.nvert;i++)
	{
	var vi=this.vertices[i];
	for (var j=0;j<this.dim;j++) c[j]+=vi[j]/this.nvert;
	};
		
return c;	
};

//	scale()
//	=======

//	Rescale the polygon, keeping its centre of mass constant

GonND.prototype.scale =
function (factor)
{
var c=com();
		
//	Rescale
		
for (var i=0;i<this.nvert;i++)
	{
	var vi=this.vertices[i];
	for (var j=0;j<this.dim;j++) vi[j]=c[j]+factor*(vi[j]-c[j]);
	};
};

//	convex()
//	========

//	Check whether polygon is convex or not.

//	If the inner product of the wedge products <(e1^e2),(e2^e3)> is positive,
//	for all consecutive triples of edges, polygon is convex.
//	For 3 dimensions, that's just a dot product of cross products, (e1 x e2).(e2 x e3)

GonND.prototype.convex -
function ()
{
if (this.convexity<0)
	{
	var con=true;

	if (this.nvert>3)
		{
		var edges=this.edgeVectors();
		for (var i=0;i<this.nvert;i++)
			{
			var e1=edges[i], e2=edges[(i+1)%this.nvert], e3=edges[(i+2)%this.nvert];
			var iwp=0.0;
			for (var j=0;j<this.dim;j++)
				for (var k=j+1;k<this.dim;k++)
					iwp+=(e1[j]*e2[k]-e1[k]*e2[j])*(e2[j]*e3[k]-e2[k]*e3[j]);

			if (iwp<0)
				{
				con=false;
				break;
				};
			};
		};
	this.convexity=con?1:0;
	};
return (this.convexity==1);
};

//	edgeVectors()
//	=============

//	Returns an array of edge vectors

GonND.prototype.edgeVectors =
function ()
{
var ev=Utils.multiDim([this.nvert,this.dim],0);
for (var i=0;i<this.nvert;i++)
	{
	var v0=this.vertices[i], v1=this.vertices[(i+1)%this.nvert], evi=ev[i];
	for (var k=0;k<this.dim;k++) evi[k]=v1[k]-v0[k];
	};
return ev;
};

//	setEdgeIDs()
//	============

//	Set all the edge IDs to current edge numbers, with EDGE_ORIGINAL flag

GonND.prototype.setEdgeIDs =
function ()
{
this.edgeIDs=new Array(this.nvert);
for (var i=0;i<this.nvert;i++) this.edgeIDs[i]=GonND.EDGE_ORIGINAL | i;
};

//	splitLinks()
//	============

//	Identify the line segment(s) that comprise the intersection of this GonND
//	with a line/plane/hyperplane (the SPLITTING SET), to a specified tolerance.
//
//	"svec[0 ... dim-1]" give a unit vector normal to the splitting set;
//	"svec[dim]" gives the distance of the splitting set from the origin in the
//	direction of the normal (i.e. it is the dot product of the normal with any point in
//	the splitting set).
//
//	"tolerance" specifies how far from the splitting set a point can be, for it to be
//	treated as lying exactly on it.
//
//	If "pm" is not null, it gives a PointManager used for the creation of any new
//	points, avoiding duplication.
//
//	sig is a caller-supplied array of size nvert, which is set to -1/0/+1 for
//	each vertex, indicating which side of the splitting set it lies.
//
//	eLcode and vLcode are caller-supplied arrays of size nvert, which are
//	set to codes indicating the *other* point to which edge and vertex points on
//	the boundary are linked.
//	These codes are:	1+edge no.			for point on split edge
//						-1-vertex no.		for boundary vertex
//						0					no link
//
//	spts is a caller-supplied array of size nvert, which is set to the split
//	points created on the edges, or left as null when the edge is unsplit.
//
//	Value returned:
//
//		-1	whole of object lies on the -ve side
//		 0	whole of object lies in the splitting set
//		+1	whole of object lies on the +ve side
//		 2	object crosses the splitting set

GonND.prototype.splitLinks =
function (svec,	tolerance, pm, sig, eLcode, vLcode, spts)
{
//	Compute distances from splitting set for all vertices.

var sd=-svec[this.dim];
var dP=0, dN=0, dZ=0;				//	Counts of +ve, -ve, (roughly) zero distances
var dist=new Array(this.nvert);
for (var i=0;i<this.nvert;i++)
	{
	var d=sd, vi=this.vertices[i];
	for (var k=0;k<this.dim;k++) d+=svec[k]*vi[k];
	dist[i]=d;
	if (Math.abs(d)<=tolerance) {sig[i]=0; dZ++;}
	else if (d<0) {sig[i]=-1; dN++;}
	else {sig[i]=1; dP++;}
	};
	
//	Eliminate any "saw-tooth grazes"
//
//	\  /
//   \/
//	----

if (this.nvert>3 && dZ!=0 && dZ !=this.nvert)
	{
	for (var i=0;i<this.nvert;i++)
		{
		var j=(i+1)%this.nvert, k=(i+2)%this.nvert, si;
		if (sig[j]==0 && (si=sig[i])==sig[k])
			{
			sig[j]=si;
			dZ--;
			if (si==1) dP++; else dN++;
			};
		};
	};
	
//	Handle non-splitting cases

if (dZ==this.nvert)			//	All on set, or near enough
	{
	return 0;
	}
else if (dN+dZ==this.nvert)		//	All on -ve side, or near enough
	{
	return -1;
	}
else if (dP+dZ==this.nvert)		//	All on +ve side, or near enough
	{
	return 1;
	};
	
//	Handle splitting cases

if (eLcode==null) return 2;

//	Collect all boundary points in bpts, with codes for where they came
//	from in bpCode
//
//	These codes are:	1+edge no.			for point on split edge
//						-1-vertex no.		for boundary vertex

var bpts=new Array(2*this.nvert);
var bpCode=new Array(2*this.nvert), nbp=0;

//	Create points on edges that are split

var sigi;
for (var i=0;i<this.nvert;i++)
	{
	if ((sigi=sig[i])==0)
		{
		bpts[nbp]=this.vertices[i];
		bpCode[nbp++]=-1-i;
		continue;
		};

	var j=(i+1)%this.nvert;
	if (sigi*sig[j]==-1)	//	Edge crosses boundary
		{
		//	Create new point on edge, where it crosses boundary

		var di=dist[i], lambda=di/(di-dist[j]), oml=1.0-lambda;
		var newpt=new Array(this.dim);
		var vi=this.vertices[i], vj=this.vertices[j];
		for (var k=0;k<this.dim;k++) newpt[k]=lambda*vj[k]+oml*vi[k];
		if (!pm) spts[i]=newpt;
		else spts[i]=pm.findOrAdd(newpt,0).pt;
		bpts[nbp]=spts[i];
		bpCode[nbp++]=1+i;
		};
	};
	
if (nbp<2) return (dP>0 ? 1: -1);
	
//	Parameterise boundary points

var bpar=new Array(nbp);
var b0=bpts[0], b1=bpts[1], bvec=new Array(this.dim);
for (var k=0;k<this.dim;k++) bvec[k]=b1[k]-b0[k];

for (var i=0;i<nbp;i++) bpar[i]=Utils.dot(bvec,bpts[i]);
	
//	Sort boundary point codes by parameter;
//	(no need to sort pointers to coordinates)

GonND.boundarySort(bpar, bpCode, nbp);
	
//	Boundary points are linked in pairs, with the link codes in
//	eLcode and vLcode for boundary points on split edges,
//	and vertices that lie on the boundary.

var crossIn=-1;

var inGon=false;
for (var i=0;i<nbp;i++)
	{
	var bpCi=bpCode[i];
	
	//	Normally we just invert inGon at each boundary point; however,
	//	we forbid in-crossings if the next boundary point is a vertex
	//	joined to the current one; this avoids mistaking an edge running
	//	along the boundary for entry to the interior.

	if (!inGon && i+1<nbp && bpCi<0 && bpCode[i+1]<0)
		{
		var c1=-1-bpCi, c2=-1-bpCode[i+1];
		if (c1==(c2+1)%this.nvert || c2==(c1+1)%this.nvert) continue;
		};
	
	inGon=!inGon;
	if (inGon) crossIn=bpCi;	//	Remember crossing in, to use when we cross out
	else						//	Mutually link crossing-in and crossing-out points
		{
		if (crossIn<0) vLcode[-1-crossIn]=bpCi; else eLcode[crossIn-1]=bpCi;
		if (bpCi<0) vLcode[-1-bpCi]=crossIn;	else eLcode[bpCi-1]=crossIn;
		};
	};
	
return 2;
}

//	splitSegs()
//	============

//	Returns the line segment(s) that comprise the intersection of this GonND
//	with a line/plane/hyperplane (the SPLITTING SET), to a specified tolerance.
//
//	"svec[0 ... dim-1]" give a unit vector normal to the splitting set;
//	"svec[dim]" gives the distance of the splitting set from the origin in the
//	direction of the normal (i.e. it is the dot product of the normal with any point in
//	the splitting set).
//
//	"tolerance" specifies how far from the splitting set a point can be, for it to be
//	treated as lying exactly on it.
//
//	If "pm" is not null, it gives a PointManager used for the creation of any new
//	points, avoiding duplication.
//
//	If "includeEdges" is true, edges of the GonND that graze the splitting set are
//	included.
//
//	"isEdge"[0] is set equal to an array of flags identifying grazing edges, if
//	it's not null, and if "includeEdges" is true
//
//	The array returned contains 0 or more pairs of points, which are the line
//	segments' endpoints.

GonND.prototype.splitSegs =
function (svec, tolerance, pm, includeEdges, isEdge)
{
var sig=Utils.multiDim([this.nvert],0), eLcode=Utils.multiDim([this.nvert],0), vLcode=Utils.multiDim([this.nvert],0);
var spts=Utils.multiDim([this.nvert],null);
var slr=this.splitLinks(svec, tolerance, pm, sig, eLcode, vLcode, spts);

//	Count up boundary points, and maybe count grazing edges

var nbp=0, nsegs=0;
for (var i=0;i<this.nvert;i++)
	{
	if (eLcode[i]!=0 || vLcode[i]!=0) nbp++;
	if (includeEdges && sig[i]==0 && sig[(i+1)%this.nvert]==0) nsegs++;
	};
nsegs+=nbp/2;

var segs=Utils.multiDim([nsegs,2],null);

if (includeEdges && isEdge!=null) isEdge[0]=new Array(nsegs);

//	Transcribe segments starting from a vertex

var iseg=0, lc;
for (var i=0;i<this.nvert;i++)
	{
	if ((lc=vLcode[i])!=0)
		{
		segs[iseg][0]=this.vertices[i];
		vLcode[i]=0;
		if (lc<0) {segs[iseg][1]=this.vertices[-1-lc]; vLcode[-1-lc]=0;}
		else {segs[iseg][1]=spts[lc-1]; eLcode[lc-1]=0;};
		iseg++;
		};
	};
	
//	Transcribe segments starting from mid-edge

for (var i=0;i<this.nvert;i++)
	{
	if ((lc=eLcode[i])!=0)
		{
		segs[iseg][0]=spts[i];
		eLcode[i]=0;
		if (lc<0) {segs[iseg][1]=this.vertices[-1-lc]; vLcode[-1-lc]=0;}
		else {segs[iseg][1]=spts[lc-1]; eLcode[lc-1]=0;};
		iseg++;
		};
	};
	
//	Transcribe grazing edges

if (includeEdges)
for (var i=0;i<this.nvert;i++)
if (sig[i]==0 && sig[(i+1)%this.nvert]==0)
	{
	segs[iseg][0]=this.vertices[i];
	segs[iseg][1]=this.vertices[(i+1)%this.nvert];
	if (isEdge!=null) isEdge[0][iseg]=true;
	iseg++;
	};

return segs;
}

//	BSPsplit()
//	==========

//	Split this BSPobject into parts that lie on either side of, or within, a specified
//	line/plane/hyperplane (the SPLITTING SET), to a specified tolerance.
//
//	"svec[0 ... dim-1]" give a unit vector normal to the splitting set;
//	"svec[dim]" gives the distance of the splitting set from the origin in the
//	direction of the normal (i.e. it is the dot product of the normal with any point in
//	the splitting set).
//
//	"ID" is an ID number associated with the choice of splitting set.
//
//	"sContents" is a Vector containing BSPobjects in the splitting set, or null.  These
//	objects might be used to further modify the parts created.
//
//	"tolerance" specifies how far from the splitting set a point can be, for it to be
//	treated as lying exactly on it.
//
//	If "pm" is not null, it gives a PointManager used for the creation of any new
//	points, avoiding duplication.
//
//	The parent BSPobject itself might be placed in one of these Vectors, if it lies
//	wholly on one side of (or wholly within) the splitting set.  Supplying null
//	for these arguments means a value is returned, but no actual splitting is done.
//
//	Value returned:
//
//		-1	whole of object lies on the -ve side
//		 0	whole of object lies in the splitting set
//		+1	whole of object lies on the +ve side
//		 2	object crosses the splitting set

//	The GonND objects created all have "surfaceID" and "flags" inherited from the parent,
//	and "edgeIDs" arrays with edge identifiers either:
//	 -	inherited, maybe with EDGE_CUT_AT_START or EDGE_CUT_AT_END set.
//	 -	set to ID, with EDGE_POSITIVE or EDGE_NEGATIVE bit set, and with EDGE_IN_CONTENT
//		set if sContents is not null, and the edge lies wholly within one of the objects
//		in SContents.
//
//	(When the parent GonND has a null edgeIDs array, so do the pieces they split to.)

GonND.prototype.BSPsplit =
function (svec, ID, sContents, tolerance, pm, negSplit, onSplit, posSplit)
{
var sig=Utils.multiDim([this.nvert],0), eLcode=Utils.multiDim([this.nvert],0), vLcode=Utils.multiDim([this.nvert],0);
var spts=Utils.multiDim([this.nvert],null);
var res=this.splitLinks(svec, tolerance, pm, sig, eLcode, vLcode, spts);
	
//	Handle non-splitting cases

if (res!=2)
	{
	if (res==0)			//	All on set, or near enough
		{
		if (onSplit!=null) onSplit.push(this);
		}
	else if (res==-1)		//	All on -ve side, or near enough
		{
		if (negSplit!=null) negSplit.push(this);
		}
	else if (res==1)		//	All on +ve side, or near enough
		{
		if (posSplit!=null) posSplit.push(this);
		};
	return res;
	};
	
//	Handle splitting cases

if (posSplit==null || negSplit==null) return 2;

//	Peel off new GonNDs that lie wholly on one side or the other

var ngv=new Array(2*this.nvert);
var haveEIDS = (this.edgeIDs != null);
var eids = haveEIDS ? new Array(2*this.nvert) : null;
var split;
var side=0;
var used=Utils.multiDim([this.nvert],false);

for (var qq=0;qq<10000;qq++)
	{
	
	//	Set p0 to index of first unused non-boundary vertex
	
	var p0;
	for (p0=0;p0<this.nvert;p0++)
		if ((side=sig[p0])!=0 && !used[p0])
			break;
	if (p0==this.nvert) break;				//	Used them all up
	
	var nv=0;							//	Count piece's vertices
	var cp=p0;							//	Current parent vertex index
	
	var eid=ID | (side<0?GonND.EDGE_NEGATIVE:GonND.EDGE_POSITIVE);
	
	for (var rr=0;rr<10000;rr++)
		{
		var lc=0;						//	Set to link code if we hit boundary
		
		ngv[nv]=this.vertices[cp];			//	Copy current parent vertex
		if (haveEIDS) eids[nv]=this.edgeIDs[cp];	//	Copy the following edge ID
		nv++;
		used[cp]=true;					//	Ensure we don't try to reuse
		
		if ((split=spts[cp])!=null)		//	Edge we're on is split
			{
			ngv[nv]=split;				//	Copy the split-point
			if (haveEIDS)
				{
				eids[nv-1]|=GonND.EDGE_CUT_AT_END;	//	Mark previous edge as cut
				eids[nv]=eid;			//	Split-created edge
				};
			nv++;
			lc=eLcode[cp];				//	See what it links to
			}
		else
			{
			cp=(cp+1)%this.nvert;			//	Next vertex point
			if (sig[cp]==0 && vLcode[cp]!=0)	//	... lies on boundary
				{
				ngv[nv]=this.vertices[cp];	//	Copy boundary vertex point
				if (haveEIDS) eids[nv]=eid;	//	Split-created edge
				nv++;
				lc=vLcode[cp];			//	See what it links to
				};
			};

		if (lc!=0)
			{
			if (lc<0)
				{
				cp=-1-lc;					//	Link to a boundary vertex point
				}
			else
				{
				cp=lc-1;					//	Link to a mid-edge point
				ngv[nv]=spts[cp];			//	Copy that point
				if (haveEIDS)
					eids[nv]=this.edgeIDs[cp]|GonND.EDGE_CUT_AT_START;
											//	Copy the following edge ID
				nv++;
				cp=(cp+1)%this.nvert;			//	Move on to next vertex
				};
			};

		if (cp==p0) break;				//	Completed loop
		};
		
	var piece=this.newGon(ngv,0,nv);
	piece.surfaceID=this.surfaceID;
	piece.flags=this.flags;
	if (haveEIDS)
		{
		piece.edgeIDs=new Array(nv);
		for (var k=0;k<nv;k++) piece.edgeIDs[k]=eids[k];
		if (sContents!=null) piece.markEdges(svec, sContents, eid, tolerance, pm);
		};
	(side<0?negSplit:posSplit).push(piece);
	};
return 2;
}

//	Sort boundary points by parameter, taking accompanying code along with them

GonND.boundarySort =
function (bpar, bpCode, nbp)
{
for (var i=0;i<nbp-1;i++)
	{
	var minpar=Number.POSITIVE_INFINITY;
	var mj=-1;
	for (var j=i;j<nbp;j++)		//	Locate minimum in set from i -> nbp-1
		{
		if (bpar[j]<minpar)
			{
			minpar=bpar[j];
			mj=j;
			};
		};
	if (mj!=i)					//	Swap to bottom of set
		{
		bpar[mj]=bpar[i];	bpar[i]=minpar;
		var tmp=bpCode[mj];	bpCode[mj]=bpCode[i];	bpCode[i]=tmp;
		};
	};

}


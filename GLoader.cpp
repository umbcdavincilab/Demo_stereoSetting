// ----------------------------------------------------------------------------
// CGLoader : the class implementation
//
// Creation : May 25  2014
//
// Copyright(C) 2012-2014  DaVinCI @ UMBC
//
// ----------------------------------------------------------------------------
#include "GLoader.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//  the class implementation of CTgdataLoader 
CTgdataLoader::CTgdataLoader(bool bLoadColor) :
	CGLoader < vector< GLfloat > > ("tgdata"),
	m_bWithColor(bLoadColor)
{
}

CTgdataLoader::~CTgdataLoader()
{
}

int CTgdataLoader::load(const string& fn)
{
	int ret = CGLoader< vector< GLfloat> >::load(fn);
	if ( 0 != ret ) {
		return ret;
	}

	unsigned long nLine = 0, icntLine = 0;
	unsigned long nPoints = 0, icntPoint = 0, nPtTotal = 0;
	int nStride = m_bWithColor? 6:3;

	GLfloat d;
	unsigned long icurFline = 0;
	istringstream istr(m_pdata, istringstream::in);
	istr >> nLine;
	icurFline ++;

	if ( nLine < 1 ) {
		cerr << "empty source data." << endl;
		return -3;
	}

	try {
		m_pbuf = new vector<GLfloat> [ nLine ];
	}
	catch (std::bad_alloc & e) {
		cerr << "FATAL: alloc failed for the head list - " << e.what() << endl;
		return -4;
	}

	/**** WE PRESUME THERE IS NO EMPTY LINE IN THE SOURCE FILE ********/

	// iterate through all polylines
	while (istr.good() && icntLine < nLine) {
		istr >> nPoints;
		// nPoints must be integer
		if ( nPoints != (unsigned long)(nPoints) || nPoints >= 0xffffffff) {
			goto err_parse;
		}
		icurFline ++;
		icntPoint = 0;
		m_pbuf [ icntLine ].resize(nStride*nPoints);
		nPtTotal += nPoints;

		// iterate through all points for a single polyline
		while (istr.good() && icntPoint < nPoints) {

			// read info for a single point once
			for (int icnt = 0, pos = 0; istr && icnt < nStride; icnt++) {
				pos = icnt + icntPoint*nStride;
				if ( m_bWithColor ) {
					pos += (icnt>2?-3:3);
				}

				istr >> fixed >> setprecision(6) >> d;
				m_pbuf [ icntLine ][ pos ] = d;
				icurFline ++;

				if ( icnt < 3 ) {
					if( m_maxCoord[ icnt ] < d ) {
						m_maxCoord[ icnt ] = d;
					}
				    if( m_minCoord[ icnt ] > d ) {
						m_minCoord[ icnt ] = d;
				    } 
				}
			}

			if ( ! istr ) {
				goto err_parse;
			}

			icntPoint ++;
		}

		if ( icntPoint < nPoints ) {
			goto err_parse;
		}

		icntLine ++;
	}

	if ( icntLine < nLine ) {
		goto err_parse;
	}

	m_numel = nLine;
	report();
	cout << nPtTotal << " vertices totally." << endl;
	unmap();
	return 0;

err_parse:
	cerr << "FATAL: parsing source failed at Line No." << icurFline << 
		" -- File not compatible with designated format: " << 
		g_gtypes[m_gtype] << endl;
	unmap();
	m_numel = 0;
	delete [] m_pbuf;
	m_pbuf = NULL;
	return -1;
}

int CTgdataLoader::dump(const string& fn, const edge_flag_t* pedgeflags)
{
	ofstream ofs(fn.c_str());
	if ( ! ofs.is_open() ) {
		cerr << "failed to create file : " << fn << " for serializing." << endl;
		return -1;
	}

	ofs << fixed << setprecision(6);

	unsigned long szTotal = getSize(), szLine;
	unsigned long nStride = m_bWithColor? 6:3;
	unsigned long offset = m_bWithColor? 3:0;

	unsigned long szToDump = 0;
	if ( pedgeflags ) {
		for (size_t i=0; i<pedgeflags->size(); ++i) {
			szToDump += ( GL_TRUE == (*pedgeflags)[i] ? 1:0 );
		}
	}

	ofs << (pedgeflags?szToDump:szTotal)
		/*
		<< " lines, with" << (m_bWithColor?"":"out")
		<< " colors"	<< endl 
		*/
		<< endl;

	for (unsigned long idx = 0; idx < szTotal; ++idx) {
		// apply edge Flags (here it is streamline flag in fact) as filters
		if ( pedgeflags && GL_FALSE == (*pedgeflags)[idx] ) {
			continue;
		}
		vector<GLfloat> & curLine = getElement( idx );
		szLine = static_cast<unsigned long> ( curLine.size() );
		
		ofs << 
			/*
			"curline No." << idx << " : [" << szLine << 
			" elements, " << 
			*/
			( szLine/nStride )	//<< " points ]" 
			<< endl;
		
		for (unsigned long j = 0; j < szLine; j += nStride) {
			ofs << curLine[j+offset+0] <<  " " << 
				curLine[j+offset+1] << " " << curLine[j+offset+2];
			if ( m_bWithColor ) {
				ofs << " " << curLine[j] <<  " " << 
					curLine[j+1] << " " << curLine[j+2];
			}
			ofs << endl;
		}
	}

	ofs.close();

	return 0;
}

/////////////////////////////////////////////////////////////////////////////
//  the class implementation of CStreamtubes
CStreamtubes::CStreamtubes() :
    m_loader(),
	m_lod(5),
	m_fAdd(0.5),
	m_fRadius(0.25),
	m_fbdRadius(20.0)
{
	for (int i=0;i<3;++i) {
		m_maxCoord[i] = -LOCAL_MAXDOUBLE;
		m_minCoord[i] = LOCAL_MAXDOUBLE;
	}
}

CStreamtubes::~CStreamtubes()
{
}

int CStreamtubes::loadGeometry(const string& fn) 
{
	if ( 0 != m_loader.load(fn) ) {
		cout << "Loading geometry failed - loader aborted abnormally.\n";
		return -1;
	}

	m_fbdRadius = m_loader.getBoundRadius();

	// generate streamtube, each for a streamline
	unsigned long szTotal = m_loader.getSize();
	m_alltubevertices.resize( szTotal );
	m_alltubenormals.resize ( szTotal );
	m_alltubecolors.resize ( szTotal );
	m_alltubefaceIdxs.resize( szTotal );

	cout << "Generating streamtube meshes .... ";
	for (unsigned long idx = 0; idx < szTotal; ++idx) {
		buildTubefromLine( idx );
	}

	cout << " finished.\n";

	cout << "X : " << m_minCoord[0] << " - " << m_maxCoord[0] << "\n";
	cout << "Y : " << m_minCoord[1] << " - " << m_maxCoord[1] << "\n";
	cout << "Z : " << m_minCoord[2] << " - " << m_maxCoord[2] << "\n";
	return 0;
}

void CStreamtubes::getExtent(GLfloat x[2], GLfloat y[2], GLfloat z[2])
{
	x[0] = m_minCoord[0], x[1] = m_maxCoord[0];
	y[0] = m_minCoord[1], y[1] = m_maxCoord[1];
	z[0] = m_minCoord[2], z[1] = m_maxCoord[2];
}

void CStreamtubes::buildTubefromLine(unsigned long lineIdx) 
{
	const vector<GLfloat> & line = m_loader.getElement( lineIdx );
	unsigned long szPts = static_cast<unsigned long>( line.size()/6 );
	GLfloat x1,x2,y1,y2,z1,z2,dx,dy,dz;
	GLfloat rax, ray, raz; // for a vector acting as a rotating axis
	GLfloat angle, theta;
	GLfloat nx,ny,nz; // rotated normals
	GLfloat scaleFactor; // scaling factor for generating quads

	// geometry generated for constructing the streamtube
	vector<GLfloat>		tube_vertices;
	vector<GLfloat>		tube_normals;
	vector<GLfloat>		tube_colors;
	vector<GLuint>		tube_faceIdxs;

	/* vertex array for each tube segment will contain all vertices needed to
	 * contruct the faces, tinted by the color originally for the start point
	 * of the source streamline segment,which is stored in the color array
	 */
	tube_vertices.resize ( szPts * m_lod * 3 );
	tube_normals.resize( szPts * m_lod * 3 );

	// szPts points contain szPts-1 line segments
	//tube_colors.resize( (szPts-1) * 3 );
	tube_colors.resize( szPts * m_lod * 3 );
	tube_faceIdxs.resize( (szPts-1) * m_lod * 4 );

	for (unsigned long idx1 = 0, idx2 = 1; idx1 < szPts; idx1++) {
		idx2 = idx1 + 1;
		if ( szPts-1 == idx1 ) {
			idx2 = idx1 - 1;
		}

		/* for each segment in a single streamline, try to interpolate
		 * a series of extraneous points to make the line segment look like a
		 * tube segment by contructing a ring around the line fragment
		 */
		x1 = line [ idx1*6 + 3 ], 
		   y1 = line [ idx1*6 + 4 ], 
		   z1 = line [ idx1*6 + 5 ];

		x2 = line [ idx2*6 + 3 ], 
		   y2 = line [ idx2*6 + 4 ], 
		   z2 = line [ idx2*6 + 5 ];

		scaleFactor = m_fRadius;

		// map the streamline segment vector's coordinate to the RGB color
		// space
		if ( szPts-1 == idx1 ) {
			dx = x1 - x2, dy = y1 - y2, dz = z1 - z2;
		}
		else {
			dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
		}

		/* firstly project a third vector in between each pair of adjacent line
		 * fragments
		 */
		if ( 0 == idx1) {
			dx = x1 - x2, dy = y1 - y2, dz = z1 - z2;
		}
		else if ( szPts-1 == idx1 ) {
			dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
		}
		else { // idx1 >= 1
			GLfloat x0 = line [ (idx1-1)*6 + 3 ],
					y0 = line [ (idx1-1)*6 + 4 ],
					z0 = line [ (idx1-1)*6 + 5 ];
			GLfloat dx1 = x0 - x1, dy1 = y0 - y1, dz1 = z0 - z1,
					dx2 = x1 - x2, dy2 = y1 - y2, dz2 = z1 - z2;
			normalize(dx1, dy1, dz1);
			normalize(dx2, dy2, dz2);

			dx = m_fAdd*dx1 + dx2, dy = m_fAdd*dy1 + dy2, dz = m_fAdd*dz1 + dz2;
		}

		normalize(dx, dy, dz);

		/* decide the axis and angle used for rotating
		 * the angle between axis and the third vector interpolated above
		 * which is calculated as  "acos ( dotproduct( (0,0,1),(dx,dy,dz) ) )"
		 */
		angle = acos (  dz );
		/*
		 * the axis is simply the vector perpendicular both to these two vectors
		 * i.e. (rax,ray,raz) = crossproduct( (0,0,1),(dx,dy,dz) ).
		 */
		rax = -dy, ray = dx, raz = 0;

		normalize(rax, ray, raz);

		/* construct a series of unit normals in the xy plane, then rotate each of them
		 * around the axis by the angle above so as to contruct a multiple of ring that
		 * is expected to appear as a tube
		 */
		for (GLubyte l = 0; l < m_lod; ++l) {
			theta = 2 * 3.1415926 / m_lod * l;

			// rotate a unit normal ( cos(theta), sin(theta), 0 ) by the degree of
			// angle, around the axis above, then store the resulting normals
			tube_normals[ idx1*3 * m_lod + 3*l + 0 ] = nx = 
				( (1 - cos(angle)) * rax * rax + cos(angle) ) * cos(theta) +
				( (1 - cos(angle)) * rax * ray - sin(angle)*raz ) * sin(theta) ;

			tube_normals[ idx1*3 * m_lod + 3*l + 1 ] = ny =
				( (1 - cos(angle)) * rax * ray + sin(angle)*raz ) * cos(theta) +
				( (1 - cos(angle)) * ray * ray + cos(angle) ) * sin(theta);

			tube_normals[ idx1*3 * m_lod + 3*l + 2 ] = nz =
				( (1 - cos(angle)) * rax * raz - sin(angle)*ray ) * cos(theta) +
				( (1 - cos(angle)) * ray * raz + sin(angle)*rax ) * sin(theta) ;

			tube_vertices[ idx1*3 * m_lod + 3*l + 0 ] = scaleFactor*nx + line[idx1*6+3];
			tube_vertices[ idx1*3 * m_lod + 3*l + 1 ] = scaleFactor*ny + line[idx1*6+4];
			tube_vertices[ idx1*3 * m_lod + 3*l + 2 ] = scaleFactor*nz + line[idx1*6+5];
			tube_colors [ idx1*3 * m_lod + 3*l + 0 ] = line [ idx1*6 + 0 ];
			tube_colors [ idx1*3 * m_lod + 3*l + 1 ] = line [ idx1*6 + 1 ];
			tube_colors [ idx1*3 * m_lod + 3*l + 2 ] = line [ idx1*6 + 2 ];

			// find the maximal and minimal coordinats among the new
			// vertices of the streamtubes
			for (int j=0; j<3; ++j) {
				if ( tube_vertices[ idx1*3 * m_lod + 3*l + j ] > m_maxCoord[j] ) {
					m_maxCoord[j] = tube_vertices[ idx1*3 * m_lod + 3*l + j ];
				}
				if ( tube_vertices[ idx1*3 * m_lod + 3*l + j ] < m_minCoord[j] ) {
					m_minCoord[j] = tube_vertices[ idx1*3 * m_lod + 3*l + j ];
				}
			}

			// the tube is finally established by a multiple of quads
			// and here we use vertex index to represent faces, each for a quad
			if ( szPts-1 > idx1 ) {
				tube_faceIdxs [ idx1*4 * m_lod + 4*l + 0 ] = idx1*m_lod + l;
				tube_faceIdxs [ idx1*4 * m_lod + 4*l + 1 ] = (idx1+1)*m_lod + l;
				tube_faceIdxs [ idx1*4 * m_lod + 4*l + 2 ] = (idx1+1)*m_lod + (l + 1)%m_lod;
				tube_faceIdxs [ idx1*4 * m_lod + 4*l + 3 ] = idx1*m_lod + (l + 1)%m_lod;
			}
		}
	}

	// save the streamtube geometry for current streamline
	m_alltubevertices[ lineIdx ] = tube_vertices;
	m_alltubenormals[ lineIdx ] = tube_normals;
	m_alltubecolors[ lineIdx ] = tube_colors;
	m_alltubefaceIdxs[ lineIdx ] = tube_faceIdxs;
}

void CStreamtubes::draw()
{
	glPushMatrix();
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_COLOR_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );

	glEnable(GL_NORMALIZE);
	glEnable(GL_RESCALE_NORMAL);
	glFrontFace(GL_CCW);
	glEnable(GL_CULL_FACE);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// by default, enable lighting and use smooth shading
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glShadeModel(GL_SMOOTH);

	// move the local coordinate system so that the default origin is still
	// located at the center of the object coordinate system
	glTranslatef( -( m_minCoord[0] + m_maxCoord[0] )/2,
			-( m_minCoord[1] + m_maxCoord[1] )/2,
			-( m_minCoord[2] + m_maxCoord[2] )/2);
	/* load streamtube geometry and render
	*/
	unsigned long szTotal = m_loader.getSize();
	for (unsigned long idx = 0; idx < szTotal; ++idx) {

		glVertexPointer(3, GL_FLOAT, 0, &m_alltubevertices[idx][0]);

		glColorPointer(3, GL_FLOAT, 0, &m_alltubecolors[idx][0]);

		glNormalPointer(GL_FLOAT, 0 ,&m_alltubenormals[idx][0]);

		glDrawElements(GL_QUADS, m_alltubefaceIdxs[idx].size(), 
				GL_UNSIGNED_INT, &m_alltubefaceIdxs[idx][0]);

	}
	glPopAttrib();
	glPopMatrix();
}


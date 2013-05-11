/***********************************************************
*   Simple GLSL demo code for COS426 Computer Graphics     *
*                           ****                           *
*   Most is borrowed from                                  *
*   http://www.lighthouse3d.com/opengl/glsl/               *
*                           ****                           *
*   Dependencies: glew, glut, gl                           *
*                           ****                           *
*   Pulled together by Aleksey Boyko(aboyko@princeton.edu) *
************************************************************/

#include "Utilities.h"
#include <stdio.h>
#include <string.h>

float initGlew(bool verbose)
{
	float versionSupported = 0.0;

	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		if(verbose) cerr<<"Error: "<<glewGetErrorString(err)<<endl;
	}
	else
	{
		versionSupported = (float)strtod((char*)glGetString(GL_VERSION), NULL);
		if(verbose) 
		{
			cout<<"Status: Using GLEW "<<glewGetString(GLEW_VERSION)<<endl;

			cout<<"Vendor: "<<glGetString(GL_VENDOR)<<endl;
			cout<<"Renderer:"<<glGetString(GL_RENDERER)<<endl;
			cout<<"OpenGL Version: "<<glGetString(GL_VERSION)<<endl;
			cout<<"==============\nVerifying:"<<endl;

			cout<<"OpenGL 2.0: ";
			if(glewIsSupported("GL_VERSION_2_0"))
			{
				cout<<"Supported"<<endl;
			}
			else
				cout<<"Not Supported"<<endl;
			cout<<"OpenGL 2.1: ";
			if(glewIsSupported("GL_VERSION_2_1"))
			{
				cout<<"Supported"<<endl;
				versionSupported = 2.1;
			}
			else
				cout<<"Not Supported"<<endl;
			cout<<"OpenGL 3.0: ";
			if(glewIsSupported("GL_VERSION_3_0"))
			{
				cout<<"Supported"<<endl;
				versionSupported = 3.;
			}
			else
				cout<<"Not Supported"<<endl;
			cout<<"OpenGL 3.1: ";
			if(glewIsSupported("GL_VERSION_3_1"))
			{
				cout<<"Supported"<<endl;
				versionSupported = 3.1;
			}
			else
				cout<<"Not Supported"<<endl;
			cout<<"OpenGL 3.2: ";
			if(glewIsSupported("GL_VERSION_3_2"))
			{
				cout<<"Supported"<<endl;
				versionSupported = 3.2;
			}
			else
				cout<<"Not Supported"<<endl;
			cout<<"OpenGL 3.3: ";
			if(glewIsSupported("GL_VERSION_3_3"))
			{
				cout<<"Supported"<<endl;
				versionSupported = 3.3;
			}
			else
				cout<<"Not Supported"<<endl;
			cout<<"OpenGL 4.0: ";
			if(glewIsSupported("GL_VERSION_4_0"))
			{
				cout<<"Supported"<<endl;
				versionSupported = 4.;
			}
			else
				cout<<"Not Supported"<<endl;
			cout<<"OpenGL 4.1: ";
			if(glewIsSupported("GL_VERSION_4_1"))
			{
				cout<<"Supported"<<endl;
				versionSupported = 4.1;
			}
			else
				cout<<"Not Supported"<<endl;

			cout<<"==============\nChecking available shaders:"<<endl;
			cout<<"Vertex shader: "<<(GLEW_ARB_vertex_shader?"Supported":"Not Supported")<<endl;
			cout<<"Fragment shader: "<<(GLEW_ARB_fragment_shader?"Supported":"Not Supported")<<endl;
			cout<<"Tesselation shader: "<<(GLEW_ARB_tessellation_shader?"Supported":"Not Supported")<<endl;
			cout<<"Geometry shader 4: "<<(GLEW_ARB_geometry_shader4?"Supported":"Not Supported")<<endl;
		}
	}
	return versionSupported;
}

char *textFileRead(const char *fn) {
	FILE *fp;
	char *content = NULL;

	int count=0;

	if (fn != NULL) {
		fp = fopen(fn,"rt");

		if (fp != NULL) {

			fseek(fp, 0, SEEK_END);
			count = ftell(fp);
			rewind(fp);

			if (count > 0) {
				content = (char *)malloc(sizeof(char) * (count+1));
				count = fread(content,sizeof(char),count,fp);
				content[count] = '\0';
			}
			fclose(fp);
		}
	}
	return content;
}

int textFileWrite(const char *fn, const char *s) {

	FILE *fp;
	int status = 0;

	if (fn != NULL) {
		fp = fopen(fn,"w");

		if (fp != NULL) {

			if (fwrite(s,sizeof(char),strlen(s),fp) == strlen(s))
				status = 1;
			fclose(fp);
		}
	}
	return(status);
}

void printShaderInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	int st;
	glGetShaderiv(obj,GL_COMPILE_STATUS,&st);
	cout<<"=======\nShader "<<obj<<" status: "<< (st==GL_TRUE?"OK":"Failed.")<<endl;

	glGetShaderiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetShaderInfoLog(obj, infologLength, &charsWritten, infoLog);
		printf("%s\n",infoLog);
		free(infoLog);
	}
	cout<<"======="<<endl;
}

void printProgramInfoLog(GLuint obj)
{
	int infologLength = 0;
	int charsWritten  = 0;
	char *infoLog;

	int st;
	glGetProgramiv(obj,GL_LINK_STATUS,&st);
	cout<<"=======\nProgram "<<obj<<" status: "<< (st==GL_TRUE?"OK":"Failed.")<<endl;

	glGetProgramiv(obj, GL_INFO_LOG_LENGTH,&infologLength);

	if (infologLength > 0)
	{
		infoLog = (char *)malloc(infologLength);
		glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
		printf("%s\n",infoLog);
		free(infoLog);
	}
	cout<<"======="<<endl;
}
// ***************************************************************************
// PosixThreads.c - provides a wrapper for WIN32 threading routines so that
//                  POSIX threads calls can be used.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "PosixThreads.h"

int pthread_create(pthread_t *thread, const pthread_attr_t *attr, void *(*startup)(void *), void *params) {
	DWORD threadid;
	HANDLE h;

	if(attr) h = CreateThread(attr->threadAttributes, attr->stackSize, (DWORD (WINAPI *)(LPVOID))startup, params, attr->creationFlags, &threadid);
	else h = CreateThread(NULL, 0, (DWORD (WINAPI *)(LPVOID))startup, params, 0, &threadid);

	thread->tid = threadid;

	if(!h) return -1;

	if(attr && (attr->detachState == PTHREAD_CREATE_DETACHED)) CloseHandle(h);
	else thread->handle = h;

	return 0;
}

int pthread_join(pthread_t thread, void **value_ptr) {
	DWORD ret;

	//not joinable; detached thread
	if(!thread.handle) return -1;

	ret = WaitForSingleObject(thread.handle, INFINITE);
	if(ret == WAIT_FAILED) return -1;
	else if((ret == WAIT_ABANDONED) || (ret == WAIT_OBJECT_0)) {
		if(value_ptr) GetExitCodeThread(thread.handle, (LPDWORD)value_ptr);
	}

	return 0;
}

void pthread_exit(void *value_ptr) {
	if(value_ptr)ExitThread(*(DWORD *)value_ptr);
	else ExitThread(0);
}

int pthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t *attr) {

	if(mutex) {
		if(mutex->init && !mutex->destroyed) return EBUSY;

		mutex->mutex = CreateMutex(NULL, FALSE, NULL);
		mutex->destroyed = 0;
		mutex->init = 1;
		mutex->lockedOrReferenced = 0;
	}

	return 0;
}

int pthread_mutex_lock(pthread_mutex_t *mutex) {
	DWORD ret;

	if(!mutex) return EINVAL;

	ret = WaitForSingleObject(mutex->mutex, INFINITE);

	if(ret != WAIT_FAILED) {
		mutex->lockedOrReferenced = 1;
		return 0;
	} else return EINVAL;
}

int pthread_mutex_unlock(pthread_mutex_t *mutex) {
	DWORD ret;

	if(!mutex) return EINVAL;

	ret = ReleaseMutex(mutex->mutex);

	if(ret != 0) {
		mutex->lockedOrReferenced = 0;
		return 0;
	} else return EPERM;
}

int pthread_mutex_destroy(pthread_mutex_t *mutex) {
	if(!mutex) return EINVAL;

	if(mutex->lockedOrReferenced) return EBUSY;

	mutex->destroyed = 1;
	return 0;
}

int pthread_mutexattr_init(pthread_mutexattr_t *attr) {

	if(attr) {
		attr->protocol = PTHREAD_PRIO_NONE;
		attr->pShared = PTHREAD_PROCESS_PRIVATE;
		//		attr->prioCeiling = SCHED_FIFO;
		attr->type = PTHREAD_MUTEX_DEFAULT;
		return 0;
	}

	return -1;
}

int pthread_mutexattr_destroy(pthread_mutexattr_t *attr) {
	return 0;
}

int pthread_attr_init(pthread_attr_t *attr) {
	attr->stackSize = 0;
	attr->stackAddr = NULL;
	attr->creationFlags = 0;  /*alternative is CREATE_SUSPENDED*/
	attr->threadAttributes = NULL;  /*can childthreads inherit attributes*/
	attr->detachState = PTHREAD_CREATE_JOINABLE;
	attr->contentionScope = PTHREAD_SCOPE_PROCESS;
	return 0;
}

int pthread_attr_destroy(pthread_attr_t *attr) {
	return 0;
}

int pthread_attr_setdetachstate(pthread_attr_t *attr, int detachstate) {

	/*check the validity of detach state*/
	if((detachstate != PTHREAD_CREATE_JOINABLE) && (detachstate != PTHREAD_CREATE_DETACHED))
		return EINVAL;

	/*set the deatch state*/
	if(attr) attr->detachState = detachstate;

	return 0;
}

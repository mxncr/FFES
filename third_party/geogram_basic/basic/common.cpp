/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include <geogram_basic/basic/common.h>
#include <geogram_basic/basic/process.h>
#include <geogram_basic/basic/logger.h>
#include <geogram_basic/basic/progress.h>
#include <geogram_basic/basic/command_line.h>
#include <geogram_basic/basic/stopwatch.h>
#include <geogram_basic/numerics/multi_precision.h>
#include <geogram_basic/numerics/predicates.h>
#include <sstream>
#include <iomanip>

#ifdef GEO_OS_EMSCRIPTEN
#include <emscripten.h>
#endif

namespace GEO {

    void initialize() {

        // When locale is set to non-us countries,
        // this may cause some problems when reading
        // floating-point numbers (some locale expect
        // a decimal ',' instead of a '.').
        // This restores the default behavior for
        // reading floating-point numbers.
#ifdef GEO_OS_UNIX
        setenv("LC_NUMERIC","POSIX",1);
#endif
        
        Environment* env = Environment::instance();
        // env->set_value("version", VORPALINE_VERSION);
        // env->set_value("release_date", VORPALINE_BUILD_DATE);
        // env->set_value("SVN revision", VORPALINE_SVN_REVISION);        

        Logger::initialize();
        Process::initialize();
        Progress::initialize();
        CmdLine::initialize();

        atexit(GEO::terminate);

        
        // Clear last system error
        errno = 0;

#ifdef GEO_OS_EMSCRIPTEN
        
        // This mounts the local file system when an emscripten-compiled
        // program runs in node.js.
        // Current working directory is mounted in /working,
        // and root directory is mounted in /root
        
        EM_ASM(
            if(typeof module !== 'undefined' && this.module !== module) {
                FS.mkdir('/working');
                FS.mkdir('/root');            
                FS.mount(NODEFS, { root: '.' }, '/working');
                FS.mount(NODEFS, { root: '/' }, '/root');
            }
        );
#endif        
        
    }

    void terminate() {
        if(
            CmdLine::arg_is_declared("sys:stats") &&
            CmdLine::get_arg_bool("sys:stats") 
        ) {
            Logger::div("System Statistics");
            Process::show_stats();
        }

        Progress::terminate();
        Process::terminate();
        CmdLine::terminate();
        Logger::terminate();
        Environment::terminate();
    }
}


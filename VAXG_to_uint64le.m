






  


  
  
  
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html lang="en">
  <head>
      <meta http-equiv="X-UA-Compatible" content="IE=8"/>
      
      <!-- START OF GLOBAL NAV -->
  <link rel="stylesheet" href="/matlabcentral/css/sitewide.css" type="text/css">
  <link rel="stylesheet" href="/matlabcentral/css/mlc.css" type="text/css">
  <!--[if lt IE 7]>
  <link href="/matlabcentral/css/ie6down.css" type="text/css" rel="stylesheet">
  <![endif]-->

      
      <meta http-equiv="content-type" content="text/html; charset=UTF-8">
<meta name="keywords" content="file exchange, matlab answers, newsgroup access, link exchange, matlab blog, matlab central, simulink blog, matlab community, matlab and simulink community">
<meta name="description" content="File exchange, MATLAB Answers, newsgroup access, Links, and Blogs for the MATLAB &amp; Simulink user community">
<link rel="stylesheet" type="text/css" media="print" href="/matlabcentral/css/print.css" />
<title> File Exchange - MATLAB Central</title>
<script src="/matlabcentral/fileexchange/javascripts/jquery-1.9.1.min.js?1396618518" type="text/javascript"></script>
<script src="/matlabcentral/fileexchange/javascripts/jquery-ui-1.10.1.min.js?1396618518" type="text/javascript"></script>
<script src="/matlabcentral/fileexchange/javascripts/jquery_ujs.js?1396618518" type="text/javascript"></script>
<script src="/matlabcentral/fileexchange/javascripts/application.js?1396618518" type="text/javascript"></script>
<link href="/matlabcentral/fileexchange/stylesheets/application.css?1396618518" media="screen" rel="stylesheet" type="text/css" />
<link href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later/content/WriteVAX/VAXG_to_uint64le.m" rel="canonical" />

<link rel="search" type="application/opensearchdescription+xml" title="Search File Exchange" href="/matlabcentral/fileexchange/search.xml" />


  </head>
    <body>
      <div id="header">
  <div class="wrapper">
  <!--put nothing in left div - only 11px wide shadow --> 
    <div class="main">
        	<div id="logo"><a href="/matlabcentral/?s_tid=gn_mlc_logo" title="MATLAB Central Home"><img src="/matlabcentral/images/mlclogo-whitebgd.gif" alt="MATLAB Central" /></a></div>
      
        <div id="headertools">
        

<script language="JavaScript1.3" type="text/javascript">

function submitForm(query){

	choice = document.forms['searchForm'].elements['search_submit'].value;
	
	if (choice == "entire1" || choice == "contest" || choice == "matlabcentral" || choice == "blogs"){
	
	   var newElem = document.createElement("input");
	   newElem.type = "hidden";
	   newElem.name = "q";
	   newElem.value = query.value;
	   document.forms['searchForm'].appendChild(newElem);
	      
	   submit_action = '/searchresults/';
	}
	
	switch(choice){
	   case "matlabcentral":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "c[]";
	      newElem.value = "matlabcentral";
	      document.forms['searchForm'].appendChild(newElem);
	
	      selected_index = 0;
	      break
	   case "fileexchange":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "term";
	      newElem.value = query.value;
	      newElem.classname = "formelem";
	      document.forms['searchForm'].appendChild(newElem);
	
	      submit_action = "/matlabcentral/fileexchange/";
	      selected_index = 1;
	      break
	   case "answers":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "term";
	      newElem.value = query.value;
	      newElem.classname = "formelem";
	      document.forms['searchForm'].appendChild(newElem);
	
	      submit_action = "/matlabcentral/answers/";
	      selected_index = 2;
	      break
	   case "cssm":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "search_string";
	      newElem.value = query.value;
	      newElem.classname = "formelem";
	      document.forms['searchForm'].appendChild(newElem);
	
		  submit_action = "/matlabcentral/newsreader/search_results";
	      selected_index = 3;
	      break
	   case "linkexchange":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "term";
	      newElem.value = query.value;
	      newElem.classname = "formelem";
	      document.forms['searchForm'].appendChild(newElem);
	
	      submit_action = "/matlabcentral/linkexchange/";
	      selected_index = 4;
	      break
	   case "blogs":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "c[]";
	      newElem.value = "blogs";
	      document.forms['searchForm'].appendChild(newElem);
	
	      selected_index = 5;
	      break
	   case "trendy":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "search";
	      newElem.value = query.value;
	      newElem.classname = "formelem";
	      document.forms['searchForm'].appendChild(newElem);
	
	      submit_action = "/matlabcentral/trendy";
	      selected_index = 6;
	      break
	   case "cody":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "term";
	      newElem.value = query.value;
	      newElem.classname = "formelem";
	      document.forms['searchForm'].appendChild(newElem);
	
	      submit_action = "/matlabcentral/cody/";
	      selected_index = 7;
	      break
	   case "contest":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "c[]";
	      newElem.value = "contest";
	      document.forms['searchForm'].appendChild(newElem);
	
	      selected_index = 8;
	      break
	   case "entire1":
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "c[]";
	      newElem.value = "entiresite";
	      document.forms['searchForm'].appendChild(newElem);
	      
	      selected_index = 9;
	      break
	   default:
	      var newElem = document.createElement("input");
	      newElem.type = "hidden";
	      newElem.name = "c[]";
	      newElem.value = "entiresite";
	      document.forms['searchForm'].appendChild(newElem);
	   
	      selected_index = 9;
	      break
	}

	document.forms['searchForm'].elements['search_submit'].selectedIndex = selected_index;
	document.forms['searchForm'].elements['query'].value = query.value;
	document.forms['searchForm'].action = submit_action;
}

</SCRIPT>


  <form name="searchForm" method="GET" action="" style="margin:0px; margin-top:5px; font-size:90%" onSubmit="submitForm(query)">
          <label for="search">Search: </label>
        <select name="search_submit" style="font-size:9px ">
         	 <option value = "matlabcentral">MATLAB Central</option>
          	<option value = "fileexchange" selected>&nbsp;&nbsp;&nbsp;File Exchange</option>
          	<option value = "answers">&nbsp;&nbsp;&nbsp;Answers</option>
            <option value = "cssm">&nbsp;&nbsp;&nbsp;Newsgroup</option>
          	<option value = "linkexchange">&nbsp;&nbsp;&nbsp;Link Exchange</option>
          	<option value = "blogs">&nbsp;&nbsp;&nbsp;Blogs</option>
          	<option value = "trendy">&nbsp;&nbsp;&nbsp;Trendy</option>
          	<option value = "cody">&nbsp;&nbsp;&nbsp;Cody</option>
          	<option value = "contest">&nbsp;&nbsp;&nbsp;Contest</option>
          <option value = "entire1">MathWorks.com</option>
        </select>
<input type="text" name="query" size="10" class="formelem" value="">
<input type="submit" value="Go" class="formelem gobutton" >
</form>

		  <ol id="access2">
  <li class="first">
    <a href="https://www.mathworks.com/accesslogin/createProfile.do?uri=http%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Ffileexchange%2F24793-write-vaxd-and-vaxg-files-in-r2008b-and-later%2Fcontent%2FWriteVAX%2FVAXG_to_uint64le" id="create_account_link">Create Account</a>
  </li>
  <li>
    <a href="https://www.mathworks.com/accesslogin/index_fe.do?uri=http%3A%2F%2Fwww.mathworks.com%2Fmatlabcentral%2Ffileexchange%2F24793-write-vaxd-and-vaxg-files-in-r2008b-and-later%2Fcontent%2FWriteVAX%2FVAXG_to_uint64le" id="login_link">Log In</a>
  </li>
</ol>


      </div>
	  
        <div id="globalnav">
        <!-- from includes/global_nav.html -->
        <ol>
                <li class=";" >
                        <a href="/matlabcentral/fileexchange/?s_tid=gn_mlc_fx">File Exchange</a> 
                </li>
                <li class=";" >
                        <a href="/matlabcentral/answers/?s_tid=gn_mlc_an">Answers</a> 
                </li>
                <li class=";" >
                        <a href="/matlabcentral/newsreader/?s_tid=gn_mlc_ng">Newsgroup</a> 
                </li>
                <li class=";" >
                        <a href="/matlabcentral/linkexchange/?s_tid=gn_mlc_lx">Link Exchange</a> 
                </li>
                <li class=";" >
                        <a href="http://blogs.mathworks.com/?s_tid=gn_mlc_blg">Blogs</a> 
                </li>
                <li class=";" >
                        <a href="/matlabcentral/trendy/?s_tid=gn_mlc_tnd">Trendy</a> 
                </li>
                <li class=";" >
                        <a href="/matlabcentral/cody/?s_tid=gn_mlc_cody">Cody</a> 
                </li>
                <li class=";" >
                        <a href="/matlabcentral/contest/?s_tid=gn_mlc_cn">Contest</a> 
                </li>
                <li class="icon mathworks" >
                        <a href="/?s_tid=gn_mlc_main">MathWorks.com</a> 
                </li>
        </ol>
      </div>
    </div>
  </div>
</div>

      <div id="middle">
  <div class="wrapper">

    <div id="mainbody" class="columns2">
  
  

  <div class="manifest">

    <div class="ctaBtn ctaBlueBtn btnSmall">
            <div id="download_submission_button" class="btnCont">
              <div class="btn download"><a href="/matlabcentral/fileexchange/downloads/146627/download" title="Download Now">Download Submission</a></div>
            </div>
          </div>


      <p class="license">
      Code covered by the <a href="/matlabcentral/fileexchange/view_license?file_info_id=24793" popup="new_window height=500,width=640,scrollbars=yes">BSD License</a>
      <a href="/matlabcentral/fileexchange/help_license#bsd" class="info notext preview_help" onclick="window.open(this.href,'small','toolbar=no,resizable=yes,status=yes,menu=no,scrollbars=yes,width=600,height=550');return false;">&nbsp;</a>
  </p>



  
  <h3 class="highlights_title">Highlights from <br/>
    <a href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later" class="manifest_title">Write VAXD and VAXG files in R2008b and later</a>
  </h3>
  <ul class='manifest'>
      <li class='manifest'><a href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later/content/WriteVAX/VAXD_to_uint64le.m" class="function" title="Function">VAXD_to_uint64le.m</a></li>
      <li class='manifest'><a href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later/content/WriteVAX/VAXF_to_uint32le.m" class="function" title="Function">VAXF_to_uint32le.m</a></li>
      <li class='manifest'><a href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later/content/WriteVAX/VAXG_to_uint64le.m" class="function" title="Function">VAXG_to_uint64le.m</a></li>
      <li class='manifest'><a href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later/content/WriteVAX/fwriteVAX.m" class="function" title="Function">fwriteVAX.m</a></li>
      <li class='manifest'><a href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later/content/WriteVAX/fwriteVAXD.m" class="function" title="Function">fwriteVAXD.m</a></li>
      <li class='manifest'><a href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later/content/WriteVAX/fwriteVAXG.m" class="function" title="Function">fwriteVAXG.m</a></li>
    <li class="manifest_allfiles">
      <a href="http://www.mathworks.com/matlabcentral/fileexchange/24793-write-vaxd-and-vaxg-files-in-r2008b-and-later/content/WriteVAX.zip" id="view_all_files">View all files</a>
    </li>
  </ul>

</div>


  <table cellpadding="0" cellspacing="0" class="details file contents">
    <tr>
      <th class="maininfo">
        


<div id="details">
  <h1 itemprop="name">Write VAXD and VAXG files in R2008b and later</h1>
  <p id="author">
    by 
    <span itemprop="author" itemscope itemtype="http://schema.org/Person">
      <span itemprop="name"><a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/64078">Aditya</a></span>
    </span>
  </p>
  <p>&nbsp;</p>
  <p>
    <span id="submissiondate" 
>
      20 Jul 2009
    </span>
      <span id="date_updated">(Updated 
        <span itemprop="datePublished" content="2009-07-31">31 Jul 2009</span>)
      </span>
  </p>

  <p id="summary">These functions mimic the deprecated functionality of FWRITE for writing data in VAX file formats. </p>


  
</div>

        </div>
      </th>
    </tr>
    <tr>
      <td class="file">
        <table cellpadding="0" cellspacing="0" border="0" class="fileview section">
          <tr class="title">
            <th><span class="heading">VAXG_to_uint64le.m</span></th>
          </tr>
          <tr>
            <td>
              
              <div class="codecontainer"><pre class="matlab-code">function [ uint32le ] = VAXG_to_uint64le(doubleVAXG)
%VAXG_TO_UINT64LE Converts from VAXG (double) to IEEE-LE (UINT32)
% This function converts floating point numbers initialized in MATLAB  
% into equivalent 64bit unsigned integers, and then splits the result
% into a vector columns of raw 32bit unsigned integers (little endian)
% using the specification for the VAXG floating point number format.
% VAXG is a 64bit or 'double precision' floating point format:
%  http://www.opengroup.org/onlinepubs/9629399/chap14.htm#tagfcjh_20
%
% The function is intended to be called from FWRITEVAX.
%
% See also VAXF_TO_UINT32LE, VAXD_TO_UINT64LE, FWRITEVAX
%
%  2009 The MathWorks, Inc. MATLAB and Simulink are registered trademarks
% of The MathWorks, Inc. See www.mathworks.com/trademarks for a list of 
% additional trademarks. Other product or brand names may be trademarks or 
% registered trademarks of their respective holders.

%% Define floating value properties for VAX architecture
% The generic equation for a floating point number is:
% V = (-1)^double(S) * M * A^(double(E)-B);
% Substituting M = C + F 
% V = (-1)^double(S) * (C+F) * A^(double(E)-B);
%
% Performing inverse operations to solve for E and F:
% 0 &lt;= 1 + log(M)/log(2) &lt; 1  (VAX specific)
% E = (floor) (logV / D) + 1 + B   
% F = V / ((A ^ (E-B)) - C)
%
% V = value, S = sign, M = mantissa, A = base, E = exponent, B = exponent 
% bias, C = mantissa constant, F = fraction

A = 2;      % VAX specific
B = 1024;    % VAX specific
C = 0.5;    % VAX specific
D = log(2); % VAX specific

%% Determine the sign bit. If -ve transform to positive.
S = zeros(size(doubleVAXG));
if any(doubleVAXG(:) &lt; 0)
	indices = find(doubleVAXG&lt;0);
    doubleVAXG(indices) = (-1) .* doubleVAXG(indices);
    S = zeros(size(doubleVAXG));
    S(indices) = 1;
end

%% Decompose the floating point number to SEF (Sign, Exp, Fract)
E = floor((log(doubleVAXG)./ D) + 1 + B); 
F = ((doubleVAXG ./ A.^(double(E)-B))) - C; 
% Convert floating point fraction to unsigned integer
F = F * 9007199254740992;   %VAX Specific 9007199254740992=2^53

%% Shift the bits of S, E and F
% VAX FLOAT BYTES  &lt;-----WORD1----&gt;&lt;-----WORD2----&gt;&lt;-----WORD1----&gt;&lt;-----WORD2----&gt;
% VAX FLOAT BITS   0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF
% Sign Exp Fract   SEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

S = bitshift(bitshift(uint64(S),0), 63);
E = bitshift(bitshift(uint64(E),0), 53);
F = bitshift(bitshift(uint64(F),0), 12);

%% Combine the S, E and F into the unsigned integer value
vaxInt = bitor(bitor(S,bitshift(E,-1)),bitshift(F,-12));

%% Convert to raw 32bit integers
% Flip the upper and lower bits (based on how Vax data storage format)
% VAX (64bit)      &lt;--WORD1--&gt;&lt;--WORD2--&gt;&lt;--WORD3--&gt;&lt;--WORD4--&gt;
% IEEE-LE (32bit)  &lt;--WORD2--&gt;&lt;--WORD1--&gt;
% IEEE-LE (32bit)  &lt;--WORD4--&gt;&lt;--WORD3--&gt;

% Split IEEE UINT64: &lt;--WORD1--&gt;&lt;--WORD2--&gt;|&lt;--WORD3--&gt;&lt;--WORD4--&gt;
%                    \____________________/ \____________________/
%                          vaxIntA               vaxIntB

vaxIntA = uint32(bitshift(bitshift(vaxInt,0), -32));
vaxIntB = uint32(bitshift(bitshift(vaxInt,32), -32));

% Swap WORD1 and WORD2

word2 = bitshift(bitshift(vaxIntA, 16), -16);
word1 = bitshift(vaxIntA, -16);
uint32leA = bitor(bitshift(word2,16), word1);

% Swap WORD3 and WORD4

word4 = bitshift(bitshift(vaxIntB, 16), -16);
word3 = bitshift(vaxIntB, -16);
uint32leB = bitor(bitshift(word4,16), word3);

uint32le(:,1) = reshape(uint32leA, numel(uint32leA), []);  %&lt;--WORD2--&gt;&lt;--WORD1--&gt;
uint32le(:,2) = reshape(uint32leB, numel(uint32leB), []);  %&lt;--WORD4--&gt;&lt;--WORD3--&gt;

uint32le = uint32le';

% IEEE (LE) DOUBLE BYTES  &lt;-----WORD4----&gt;&lt;-----WORD3----&gt;&lt;-----WORD2----&gt;&lt;-----WORD1----&gt;
% IEEE (LE) DOUBLE BITS   0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF
% Sign Exp Fract          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEEEEEES


</pre></div>
            </td>
          </tr>
        </table>
      </td>
    </tr>

  </table>


<p id="contactus"><a href="/company/feedback/">Contact us</a></p>

      
</div>
<div class="clearboth">&nbsp;</div>
</div>
</div>
<!-- footer.html -->
<!-- START OF FOOTER -->

<div id="mlc-footer">
  <script type="text/javascript">
function clickDynamic(obj, target_url, tracking_code) {
	var pos=target_url.indexOf("?");
	if (pos<=0) { 
		var linkComponents = target_url + tracking_code;
		obj.href=linkComponents;
	} 
}
</script>
  <div class="wrapper">
    <div>
      <ul id="matlabcentral">
        <li class="copyright first">&copy; 1994-2014 The MathWorks, Inc.</li>
        <li class="first"><a href="/help.html" title="Site Help">Site Help</a></li>
        <li><a href="/company/aboutus/policies_statements/patents.html" title="patents" rel="nofollow">Patents</a></li>
        <li><a href="/company/aboutus/policies_statements/trademarks.html" title="trademarks" rel="nofollow">Trademarks</a></li>
        <li><a href="/company/aboutus/policies_statements/" title="privacy policy" rel="nofollow">Privacy Policy</a></li>
        <li><a href="/company/aboutus/policies_statements/piracy.html" title="preventing piracy" rel="nofollow">Preventing Piracy</a></li>
        <li class="last"><a href="/matlabcentral/termsofuse.html" title="Terms of Use" rel="nofollow">Terms of Use</a></li>
        <li class="icon"><a href="/company/rss/" title="RSS" class="rssfeed" rel="nofollow"><span class="text">RSS</span></a></li>
        <li class="icon"><a href="/programs/bounce_hub_generic.html?s_tid=mlc_lkd&url=http://www.linkedin.com/company/the-mathworks_2" title="LinkedIn" class="linkedin" rel="nofollow" target="_blank"></a></li>
        <li class="icon"><a href="/programs/bounce_hub_generic.html?s_tid=mlc_fbk&url=https://plus.google.com/117177960465154322866?prsrc=3" title="Google+" class="google" rel="nofollow" target="_blank"><span class="text">Google+</span></a></li>
        <li class="icon"><a href="/programs/bounce_hub_generic.html?s_tid=mlc_fbk&url=http://www.facebook.com/MATLAB" title="Facebook" class="facebook" rel="nofollow" target="_blank"><span class="text">Facebook</span></a></li>
        		<li class="last icon"><a href="/programs/bounce_hub_generic.html?s_tid=mlc_twt&url=http://www.twitter.com/MATLAB" title="Twitter" class="twitter" rel="nofollow" target="_blank"><span class="text">Twitter</span></a></li>        
        
      </ul>
      <ul id="mathworks">
        <li class="first sectionhead">Featured MathWorks.com Topics:</li>
        <li class="first"><a href="/products/new_products/latest_features.html" onclick="clickDynamic(this, this.href, '?s_cid=MLC_new')">New Products</a></li>
        <li><a href="/support/" title="support" onclick="clickDynamic(this, this.href, '?s_cid=MLC_support')">Support</a></li>
        <li><a href="/help" title="documentation" onclick="clickDynamic(this, this.href, '?s_cid=MLC_doc')">Documentation</a></li>
        <li><a href="/services/training/" title="training" onclick="clickDynamic(this, this.href, '?s_cid=MLC_training')">Training</a></li>
        <li><a href="/company/events/webinars/" title="Webinars" onclick="clickDynamic(this, this.href, '?s_cid=MLC_webinars')">Webinars</a></li>
        <li><a href="/company/newsletters/" title="newsletters" onclick="clickDynamic(this, this.href, '?s_cid=MLC_newsletters')">Newsletters</a></li>
        <li><a href="/programs/trials/trial_request.html?prodcode=ML&s_cid=MLC_trials" title="MATLAB Trials">MATLAB Trials</a></li>
        
        		<li class="last"><a href="/company/jobs/opportunities/index_en_US.html" title="Careers" onclick="clickDynamic(this, this.href, '?s_cid=MLC_careers')">Careers</a></li>
                 
      </ul>
    </div>
  </div>
</div>
<!-- END OF FOOTER -->


      
      
<!-- SiteCatalyst code version: H.24.4.
Copyright 1996-2012 Adobe, Inc. All Rights Reserved
More info available at http://www.omniture.com -->
<script language="JavaScript" type="text/javascript" src="/scripts/omniture/s_code.js"></script>


<script language="JavaScript" type="text/javascript">



<!--
/************* DO NOT ALTER ANYTHING BELOW THIS LINE ! **************/
var s_code=s.t();if(s_code)document.write(s_code)//--></script>
<script language="JavaScript" type="text/javascript"><!--
if(navigator.appVersion.indexOf('MSIE')>=0)document.write(unescape('%3C')+'\!-'+'-')
//--></script>
<!--/DO NOT REMOVE/-->
<!-- End SiteCatalyst code version: H.24.4. -->


  

    <!-- BEGIN Qualaroo --> 
<script type="text/javascript">
  var _kiq = _kiq || [];
  (function(){
    setTimeout(function(){
    var d = document, f = d.getElementsByTagName('script')[0], s = d.createElement('script'); s.type = 'text/javascript';
    s.async = true; s.src = '//s3.amazonaws.com/ki.js/49559/ahy.js'; f.parentNode.insertBefore(s, f);
    }, 1);
  })();
</script>
<!-- END Qualaroo --> 

    </body>
</html>

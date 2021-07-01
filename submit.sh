python jemdoc.py index.jemdoc
python jemdoc.py Teaching.jemdoc
python jemdoc.py Others.jemdoc
#python jemdoc.py Research.jemdoc
python jemdoc.py Education.jemdoc
python jemdoc.py InvitedTalk.jemdoc
python jemdoc.py travel.jemdoc
python jemdoc.py Group.jemdoc
python jemdoc.py Publications_year.jemdoc
python jemdoc.py Software.jemdoc
python jemdoc.py Workshops.jemdoc
cd reading_group
python jemdoc.py index.jemdoc
python jemdoc.py paper_list.jemdoc
cd ..



scp index.html InvitedTalk.html Group.html Education.html Teaching.html Research.html travel.html Workshops.html MENU Publications_year.html Software.html jemdoc.css ../../My\ Documents/Resume/Resume_Yuxin_CHEN.pdf yc5@nobel.princeton.edu:public_html/

cd reading_group
scp index.html paper_list.html jemdoc.css yc5@nobel.princeton.edu:public_html/reading_group/

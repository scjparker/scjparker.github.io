#! /bin/sh

## edit if needed; where'e the [remote] for repository

REMOTE=origin # remote name, usually origin. Use `git remote -v` to check
CODE_BRANCH=source2 # branch where the source code will be pushed
SITE_BRANCH=master # branch where the built site will be deployed

echo "The script will build the website and deploy it to GitHub"
echo "Make sure all your changes are committed."

read -n 1 -s -r -p "Press any key to continue"

## Build: requires `bundler`
echo "\n0. Making sure dependencies are installed"
bundle install

echo "\n1. Pushing source code first"
## Push code to "site" branch in target repository
git push $REMOTE $CODE_BRANCH

echo "\n2. Building website.."
bundle exec jekyll build

echo "\n3. Pushing built website (--force)"
## git subtree deployment to $SITE_BRANCH here
## basically, push the `_site` directory to $SITE_BRANCH
git subtree split --prefix _site -b build

## Add and commit build changes
git add -A
git commit -m "Build: `date`"

git push -f $REMOTE build:$SITE_BRANCH

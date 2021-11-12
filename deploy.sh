#! /bin/sh

## edit if needed; where'e the [remote] for repository
REMOTE=origin

echo "The script will build the website and deploy it to GitHub"
echo "Make sure all your changes are committed."

read -n 1 -s -r -p "Press any key to continue"

## Build: requires `bundler`
echo "\n0. Making sure dependencies are installed"
bundle install

echo "\n1. Building website.."
bundle exec jekyll build

echo "\n2. Commiting build.."

## Add and commit build changes
git add -A
git commit -m "Build: `date`"

echo "\n3. Pushing source code"
## Push code to "site" branch in target repository
git push $REMOTE source

echo "\n4. Pushing built website (--force)"
## git subtree deployment to gh-pages here
## basically, push the `_site` directory to master
git subtree split --squash --prefix _site -b build
git push -f $REMOTE build:master
